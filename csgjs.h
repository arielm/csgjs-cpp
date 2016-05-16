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

/*
 * MODIFIED BY ARIEL MALKA (github.com/arielm) IN ORDER TO SUPPORT glm
 */

#pragma once

#include <vector>
#include <algorithm>

#include <glm/glm.hpp>

struct csgjs_vertex
{
  glm::vec3 pos;
  glm::vec3 normal;
  glm::vec2 uv;

  csgjs_vertex() = default;

  csgjs_vertex(const glm::vec3 &position, const glm::vec3 &normal, const glm::vec2 &coords)
  :
  pos(position),
  normal(normal),
  uv(coords)
  {}
};

struct csgjs_model
{
  std::vector<csgjs_vertex> vertices;
  std::vector<int> indices;
};

// public interface - not super efficient, if you use multiple CSG operations you should
// use BSP trees and convert them into model only once. Another optimization trick is
// replacing csgjs_model with your own class.

csgjs_model csgjs_union(const csgjs_model & a, const csgjs_model & b);
csgjs_model csgjs_intersection(const csgjs_model & a, const csgjs_model & b);
csgjs_model csgjs_difference(const csgjs_model & a, const csgjs_model & b);

// IMPLEMENTATION BELOW ---------------------------------------------------------------------------

// `CSG.Plane.EPSILON` is the tolerance used by `splitPolygon()` to decide if a
// point is on the plane.
static const float csgjs_EPSILON = 0.00001f;

struct csgjs_plane;
struct csgjs_polygon;
struct csgjs_node;

// Represents a plane in 3D space.
struct csgjs_plane
{
  glm::vec3 normal;
  float w;

  csgjs_plane();
  csgjs_plane(const glm::vec3 & a, const glm::vec3 & b, const glm::vec3 & c);
  bool ok() const;
  void flip();
  void splitPolygon(const csgjs_polygon & polygon, std::vector<csgjs_polygon> & coplanarFront, std::vector<csgjs_polygon> & coplanarBack, std::vector<csgjs_polygon> & front, std::vector<csgjs_polygon> & back) const;
};

// Represents a convex polygon. The vertices used to initialize a polygon must
// be coplanar and form a convex loop. They do not have to be `CSG.Vertex`
// instances but they must behave similarly (duck typing can be used for
// customization).
// 
// Each convex polygon has a `shared` property, which is shared between all
// polygons that are clones of each other or were split from the same polygon.
// This can be used to define per-polygon properties (such as surface color).
struct csgjs_polygon
{
  std::vector<csgjs_vertex> vertices;
  csgjs_plane plane;
  void flip();

  csgjs_polygon();
  csgjs_polygon(const std::vector<csgjs_vertex> & list);
};

// Holds a node in a BSP tree. A BSP tree is built from a collection of polygons
// by picking a polygon to split along. That polygon (and all other coplanar
// polygons) are added directly to that node and the other polygons are added to
// the front and/or back subtrees. This is not a leafy BSP tree since there is
// no distinction between internal and leaf nodes.
struct csgjs_csgnode
{
  std::vector<csgjs_polygon> polygons;
  csgjs_csgnode * front;
  csgjs_csgnode * back;
  csgjs_plane plane;

  csgjs_csgnode();
  csgjs_csgnode(const std::vector<csgjs_polygon> & list);
  ~csgjs_csgnode();

  csgjs_csgnode * clone() const;
  void clipTo(const csgjs_csgnode * other);
  void invert();
  void build(const std::vector<csgjs_polygon> & polygon);
  std::vector<csgjs_polygon> clipPolygons(const std::vector<csgjs_polygon> & list) const;
  std::vector<csgjs_polygon> allPolygons() const;
};
