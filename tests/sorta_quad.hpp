#ifndef SORTA_QUAD_HPP
#define SORTA_QUAD_HPP

#include "../src/geom2d.hpp"

struct sorta_quad {
  double sx, sy;
  double cx, cy;
  double ex, ey;

  sorta_quad(
    double startx = 0.0, double starty = 0.0,
    double ctrlx = 0.0, double ctrly = 0.0,
    double endx = 0.0, double endy = 0.0
  ) :sx(startx), sy(starty), cx(ctrlx), cy(ctrly), ex(endx), ey(endy)
  {}

  sorta_quad(const sorta_quad& o) :
    sx(o.sx), sy(o.sy), cx(o.cx), cy(o.cy), ex(o.ex), ey(o.ey)
  {}

  sorta_quad& rotate(double rads) {
    double c= std::cos(rads), s = std::sin(rads);
    double _sx = c*sx - s*sy, _sy = s*sx + c*sy;
    double _cx = c*cx - s*cy, _cy = s*cx + c*cy;
    double _ex = c*ex - s*ey, _ey = s*ex + c*ey;
    sx = _sx; sy = _sy;
    cx = _cx; cy = _cy;
    ex = _ex; ey = _ey;
    return *this;
  }

  sorta_quad rotated(double rads) const {
    sorta_quad ret(*this);
    ret.rotate(rads);
    return ret;
  }

  sorta_quad& translate(double tx, double ty) {
    sx+=tx; sy+=ty;
    cx+=tx; cy+=ty;
    ex+=tx; ey+=ty;
    return *this;
  }

  sorta_quad translated(double tx, double ty) const {
    sorta_quad ret(*this);
    ret.translate(tx, ty);
    return ret;
  }
};


#endif // SORTA_QUAD_HPP
