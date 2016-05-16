// Original CSG.JS library by Evan Wallace (http://madebyevan.com), under the MIT license.
// GitHub: https://github.com/evanw/csg.js/
// 
// C++ port by Tomasz Dabrowski (http://28byteslater.com), under the MIT license.
// GitHub: https://github.com/dabroz/csgjs-cpp/
// 
// Constructive Solid Geometry (CSG) is a modeling technique that uses Boolean
// operations like union and intersection to combine 3D solids. This library
// implements CSG operations on meshes elegantly and concisely using BSP trees,
// and is meant to serve as an easily understandable implementation of the
// algorithm. All edge cases involving overlapping coplanar polygons in both
// solids are correctly handled.
//

#include "csgjs.h"

using namespace std;

// Vertex implementation

// Create a new vertex between this vertex and `other` by linearly
// interpolating all properties using a parameter of `t`. Subclasses should
// override this to interpolate additional properties.
inline static csgjs_vertex interpolate(const csgjs_vertex & a, const csgjs_vertex & b, float t)
{
  csgjs_vertex ret;
  ret.pos = glm::mix(a.pos, b.pos, t);
  ret.normal = glm::mix(a.normal, b.normal, t);
  ret.uv = glm::mix(a.uv, b.uv, t);
  return ret;
}

// Plane implementation

csgjs_plane::csgjs_plane() : normal(), w(0.0f) 
{
}

bool csgjs_plane::ok() const 
{
  return glm::length(normal) > 0;
}

void csgjs_plane::flip()
{
  normal *= -1;
  w *= -1.0f;
}

csgjs_plane::csgjs_plane(const glm::vec3 & a, const glm::vec3 & b, const glm::vec3 & c)
{
  normal = glm::normalize(glm::cross(b - a, c - a));
  w = glm::dot(normal, a);
}

// Split `polygon` by this plane if needed, then put the polygon or polygon
// fragments in the appropriate lists. Coplanar polygons go into either
// `coplanarFront` or `coplanarBack` depending on their orientation with
// respect to this plane. Polygons in front or in back of this plane go into
// either `front` or `back`.
void csgjs_plane::splitPolygon(const csgjs_polygon & polygon, std::vector<csgjs_polygon> & coplanarFront, std::vector<csgjs_polygon> & coplanarBack, std::vector<csgjs_polygon> & front, std::vector<csgjs_polygon> & back) const
{
  enum
  {
    COPLANAR = 0,
    FRONT = 1,
    BACK = 2,
    SPANNING = 3
  };

  // Classify each point as well as the entire polygon into one of the above
  // four classes.
  int polygonType = 0;
  std::vector<int> types;

  for (const auto &vertex : polygon.vertices)
  {
    float t = glm::dot(normal, vertex.pos) - w;
    int type = (t < -csgjs_EPSILON) ? BACK : ((t > csgjs_EPSILON) ? FRONT : COPLANAR);
    polygonType |= type;
    types.push_back(type);
  }

  // Put the polygon in the correct list, splitting it when necessary.
  switch (polygonType) 
  {
  case COPLANAR:
    {
      if (glm::dot(normal, polygon.plane.normal) > 0)
        coplanarFront.push_back(polygon);
      else 
        coplanarBack.push_back(polygon);
      break;
    }
  case FRONT:
    {
      front.push_back(polygon);
      break;
    }
  case BACK:
    {
      back.push_back(polygon);
      break;
    }
  case SPANNING:
    {
      std::vector<csgjs_vertex> f, b;
      for (size_t i = 0; i < polygon.vertices.size(); i++) 
      {
        int j = (i + 1) % polygon.vertices.size();
        int ti = types[i], tj = types[j];
        csgjs_vertex vi = polygon.vertices[i], vj = polygon.vertices[j];
        if (ti != BACK) f.push_back(vi);
        if (ti != FRONT) b.push_back(vi);
        if ((ti | tj) == SPANNING) 
        {
          float t = (w - glm::dot(normal, vi.pos)) / glm::dot(normal, vj.pos - vi.pos);
          csgjs_vertex v = interpolate(vi, vj, t);
          f.push_back(v);
          b.push_back(v);
        }
      }
      if (f.size() >= 3) front.push_back(csgjs_polygon(f));
      if (b.size() >= 3) back.push_back(csgjs_polygon(b));
      break;
    }
  }
}

// Polygon implementation

void csgjs_polygon::flip()
{
  std::reverse(vertices.begin(), vertices.end());
  for (auto &vertex : vertices)
  {
    vertex.normal *= -1;
  }
  plane.flip();
}

csgjs_polygon::csgjs_polygon()
{
}

csgjs_polygon::csgjs_polygon(const std::vector<csgjs_vertex> & list) : vertices(list), plane(vertices[0].pos, vertices[1].pos, vertices[2].pos)
{
}

// Node implementation

// Return a new CSG solid representing space in either this solid or in the
// solid `csg`. Neither this solid nor the solid `csg` are modified.
inline static csgjs_csgnode * csg_union(const csgjs_csgnode * a1, const csgjs_csgnode * b1)
{
  csgjs_csgnode * a = a1->clone();
  csgjs_csgnode * b = b1->clone();
  a->clipTo(b);
  b->clipTo(a);
  b->invert();
  b->clipTo(a);
  b->invert();
  a->build(b->allPolygons());
  csgjs_csgnode * ret = new csgjs_csgnode(a->allPolygons());
  delete a; a = 0;
  delete b; b = 0;
  return ret;
}

// Return a new CSG solid representing space in this solid but not in the
// solid `csg`. Neither this solid nor the solid `csg` are modified.
inline static csgjs_csgnode * csg_subtract(const csgjs_csgnode * a1, const csgjs_csgnode * b1)
{
  csgjs_csgnode * a = a1->clone();
  csgjs_csgnode * b = b1->clone();
  a->invert();
  a->clipTo(b);
  b->clipTo(a);
  b->invert();
  b->clipTo(a);
  b->invert();
  a->build(b->allPolygons());
  a->invert();
  csgjs_csgnode * ret = new csgjs_csgnode(a->allPolygons());
  delete a; a = 0;
  delete b; b = 0;
  return ret;
}

// Return a new CSG solid representing space both this solid and in the
// solid `csg`. Neither this solid nor the solid `csg` are modified.
inline static csgjs_csgnode * csg_intersect(const csgjs_csgnode * a1, const csgjs_csgnode * b1)
{
  csgjs_csgnode * a = a1->clone();
  csgjs_csgnode * b = b1->clone();
  a->invert();
  b->clipTo(a);
  b->invert();
  a->clipTo(b);
  b->clipTo(a);
  a->build(b->allPolygons());
  a->invert();
  csgjs_csgnode * ret = new csgjs_csgnode(a->allPolygons());
  delete a; a = 0;
  delete b; b = 0;
  return ret;
}

// Convert solid space to empty space and empty space to solid space.
void csgjs_csgnode::invert()
{
  for (auto &polygon : polygons)
  {
    polygon.flip();
  }
  plane.flip();
  if (front) front->invert();
  if (back) back->invert();
  std::swap(front, back);
}

// Recursively remove all polygons in `polygons` that are inside this BSP
// tree.
std::vector<csgjs_polygon> csgjs_csgnode::clipPolygons(const std::vector<csgjs_polygon> & list) const
{
  if (!plane.ok()) return list;
  std::vector<csgjs_polygon> list_front, list_back;  
  for (size_t i = 0; i < list.size(); i++)
  {
    plane.splitPolygon(list[i], list_front, list_back, list_front, list_back);
  }
  if (front) list_front = front->clipPolygons(list_front);
  if (back) list_back = back->clipPolygons(list_back);
  else list_back.clear();
  
  list_front.insert(list_front.end(), list_back.begin(), list_back.end());
  return list_front;
}

// Remove all polygons in this BSP tree that are inside the other BSP tree
// `bsp`.
void csgjs_csgnode::clipTo(const csgjs_csgnode * other)
{
  polygons = other->clipPolygons(polygons);
  if (front) front->clipTo(other);
  if (back) back->clipTo(other);
}

// Return a list of all polygons in this BSP tree.
std::vector<csgjs_polygon> csgjs_csgnode::allPolygons() const
{
  std::vector<csgjs_polygon> list = polygons;
  std::vector<csgjs_polygon> list_front, list_back;
  if (front) list_front = front->allPolygons();
  if (back) list_back = back->allPolygons();
  list.insert(list.end(), list_front.begin(), list_front.end());
  list.insert(list.end(), list_back.begin(), list_back.end());
  return list;
}

csgjs_csgnode * csgjs_csgnode::clone() const
{
  csgjs_csgnode * ret = new csgjs_csgnode();
  ret->polygons = polygons;
  ret->plane = plane;
  if (front) ret->front = front->clone();
  if (back) ret->back = back->clone();
  return ret;
}

// Build a BSP tree out of `polygons`. When called on an existing tree, the
// new polygons are filtered down to the bottom of the tree and become new
// nodes there. Each set of polygons is partitioned using the first polygon
// (no heuristic is used to pick a good split).
void csgjs_csgnode::build(const std::vector<csgjs_polygon> & list)
{
  if (!list.size()) return;
  if (!plane.ok()) plane = list[0].plane;
  std::vector<csgjs_polygon> list_front, list_back;
  for (const auto &polygon : list)
  {
    plane.splitPolygon(polygon, polygons, polygons, list_front, list_back);
  }
  if (list_front.size()) 
  {
    if (!front) front = new csgjs_csgnode;
    front->build(list_front);
  }
  if (list_back.size()) 
  {
    if (!back) back = new csgjs_csgnode;
    back->build(list_back);
  }
}

csgjs_csgnode::csgjs_csgnode() : front(0), back(0)
{
}

csgjs_csgnode::csgjs_csgnode(const std::vector<csgjs_polygon> & list) : front(0), back(0)
{
  build(list);
}

csgjs_csgnode::~csgjs_csgnode()
{
  delete front;
  delete back;
}

// Public interface implementation

inline static std::vector<csgjs_polygon> csgjs_modelToPolygons(const csgjs_model & model)
{
  std::vector<csgjs_polygon> list;
  for (size_t i = 0; i < model.indices.size(); i+= 3)
  {
    std::vector<csgjs_vertex> triangle;
    for (int j = 0; j < 3; j++)
    {
      csgjs_vertex v = model.vertices[model.indices[i + j]];
      triangle.push_back(v);
    }
    list.push_back(csgjs_polygon(triangle));
  }
  return list;
}

inline static csgjs_model csgjs_modelFromPolygons(const std::vector<csgjs_polygon> & polygons)
{
  csgjs_model model;
  int p = 0;
  for (size_t i = 0; i < polygons.size(); i++)
  {
    const csgjs_polygon & poly = polygons[i];
    for (size_t j = 2; j < poly.vertices.size(); j++)
    {
      model.vertices.push_back(poly.vertices[0]);    model.indices.push_back(p++);
      model.vertices.push_back(poly.vertices[j - 1]);  model.indices.push_back(p++);
      model.vertices.push_back(poly.vertices[j]);    model.indices.push_back(p++);      
    }
  }
  return model;
}

typedef csgjs_csgnode * csg_function(const csgjs_csgnode * a1, const csgjs_csgnode * b1);

inline static csgjs_model csgjs_operation(const csgjs_model & a, const csgjs_model & b, csg_function fun)
{
  csgjs_csgnode * A = new csgjs_csgnode(csgjs_modelToPolygons(a));
  csgjs_csgnode * B = new csgjs_csgnode(csgjs_modelToPolygons(b));
  csgjs_csgnode * AB = fun(A, B);
  std::vector<csgjs_polygon> polygons = AB->allPolygons();
  delete A; A = 0;
  delete B; B = 0;
  delete AB; AB = 0;
  return csgjs_modelFromPolygons(polygons);
}

csgjs_model csgjs_union(const csgjs_model & a, const csgjs_model & b)
{
  return csgjs_operation(a, b, csg_union);
}

csgjs_model csgjs_intersection(const csgjs_model & a, const csgjs_model & b)
{
  return csgjs_operation(a, b, csg_intersect);
}

csgjs_model csgjs_difference(const csgjs_model & a, const csgjs_model & b)
{
  return csgjs_operation(a, b, csg_subtract);
}
