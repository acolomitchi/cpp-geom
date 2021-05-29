#ifndef SORTA_LINE_HPP
#define SORTA_LINE_HPP

#include <iostream>

#include "../src/geom2d.hpp"

struct sorta_line {
  double sx, sy;
  double ex, ey;
  geomalgos2d::linear_variety type;

  sorta_line(
      double x0, double y0,
      double x1, double y1,
      geomalgos2d::linear_variety t
  ) :sx(x0), sy(y0), ex(x1), ey(y1),type(t)
  {}

  bool intersection(
      const sorta_line& o,
      double& thisP, double& oP,
      double* thisIntX = nullptr, double* thisIntY = nullptr,
      double* oIntX = nullptr, double* oIntY = nullptr
  ) const
  {
    using linegeom=geomalgos2d::linegeo2d<double, double>;
    bool sameSupport = false;
    bool ret = linegeom::intersection_params(
      sx, sy, ex, ey,
      o.sx, o.sy, o.ex, o.ey,
      thisP, oP, sameSupport,
      type, o.type
    );
    using namespace geomalgos2d;
    if(ret && thisIntX)
      *thisIntX = lincomb<double>(sx, ex, thisP);
    if(ret && thisIntY)
      *thisIntY = lincomb<double>(sy, ey, thisP);
    if(ret && oIntX)
      *oIntX = lincomb<double>(o.sx, o.ex, oP);
    if(ret && oIntY)
      *oIntY = lincomb<double>(o.sy, o.ey, oP);
    return ret;
  }
};

std::ostream& operator<<(std::ostream& o, const sorta_line& l) {
    o << "[" << l.sx <<","<<l.sy<< "-" << l.type << "->"
      << l.ex << "," << l.ey << "]";
    return o;
}

#endif // SORTA_LINE_HPP
