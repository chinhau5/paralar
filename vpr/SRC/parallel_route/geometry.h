#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <boost/numeric/interval.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using point = typename bg::model::point<int, 2, bg::cs::cartesian>; 
using box = typename bg::model::box<point>;
using multi_point = typename bg::model::multi_point<point>;
using segment = typename bg::model::segment<point>;
using polygon = typename bg::model::polygon<point>;

#endif
