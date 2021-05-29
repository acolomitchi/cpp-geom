#ifndef LINEAR_GEOM_HPP
#define LINEAR_GEOM_HPP

#include <assert.h>
#include <algorithm>

#include "prec2d.hpp"
#include "dist2d.hpp"

namespace geomalgos2d {

enum linear_variety {
  line, // infinite at both sides
  ray, // semiline, start is source, no point before start, infinite at the end
  segment // limited at both sides
};

template <typename coord_t, typename al_coord_t>
    coord_t lincomb(coord_t start, coord_t end, coord_t param)
{
  al_coord_t s = static_cast<al_coord_t>(start);
  al_coord_t e = static_cast<al_coord_t>(end);
  al_coord_t p = static_cast<al_coord_t>(param);
  al_coord_t ret = (1-p)*s + p*e;
  return static_cast<coord_t>(ret);
}

template <typename coord_t>
    coord_t lincomb(coord_t start, coord_t end, coord_t param)
{
  coord_t ret = (1-param)*start + param*end;
  return ret;
}

template <typename coord_t, typename al_coord_t = double> class linegeo2d {
public:

private:

  static bool p_on_dirline(
    al_coord_t anchorX, al_coord_t anchorY,
    al_coord_t dirX, al_coord_t dirY,
    al_coord_t x, al_coord_t y,
    al_coord_t& lineP, linear_variety lineType,
    const eps_prec<al_coord_t>& prec = eps_prec<al_coord_t>::DefaultPrec
  ) {
    bool ret = false;
    bool vert = std::abs(dirX) <= prec.eps();
    bool horiz = std::abs(dirY) <= prec.eps();
    al_coord_t p=0.0;
    if (vert && horiz) { // degenerated line
      // only if the point is the same with the anchor
      ret = eps_prec<al_coord_t>::same_point_test_sq(anchorX, anchorY, x, y, prec);
      if (ret) p = 0.0;
    }
    else if (vert) {
      ret = (std::abs(x - anchorX) <= prec.eps()); // on the same vertical
      if (ret) { // x on the same vertical, safe to rely on y
        p = (y - anchorY) / dirY;
      }
    }
    else if (horiz) { // is not vert
      ret = (std::abs(y - anchorY) < prec.eps());
      if (ret) { // y on the same horizontal we can rely on x
        p = (x - anchorX) / dirX;
      }
    }
    else { // non degenerated line
      al_coord_t tx = (x - anchorX) / dirX;
      al_coord_t ty = (y - anchorY) / dirY;
      ret = (std::abs(tx - ty) <= prec.eps());
      if (ret) {
        p = (tx + ty) / 2;
      }
    }
    if(ret) {
      switch(lineType) {
        case linear_variety::segment:
          ret=(p>=0 && p<=1);
          break;
        case linear_variety::ray:
          ret = (p>=0);
          break;
        default:
          ret=true;
          break;
      }
      if(ret) lineP = p;
    }
    return ret;
  }

public:

  static bool same_point_test(
    coord_t x0, coord_t y0, coord_t x1, coord_t y1,
    const eps_prec<al_coord_t>& prec = eps_prec<al_coord_t>::DefaultPrec
  ) {
    return same_point_test_sq(x0, y0, x1, y1, prec);
  }

  static coord_t dotp(
    coord_t commonX, coord_t commonY,
    coord_t x0, coord_t y0,
    coord_t x1, coord_t y1
  ) {
    al_coord_t dx0 = static_cast<al_coord_t> (x0) - static_cast<al_coord_t> (commonX);
    al_coord_t dy0 = static_cast<al_coord_t> (y0) - static_cast<al_coord_t> (commonY);
    al_coord_t dx1 = static_cast<al_coord_t> (x1) - static_cast<al_coord_t> (commonX);
    al_coord_t dy1 = static_cast<al_coord_t> (y1) - static_cast<al_coord_t> (commonY);
    al_coord_t ret = dx0 * dx1 + dy0*dy1;
    return static_cast<coord_t> (ret);
  }

  static coord_t crossp(
    coord_t commonX, coord_t commonY,
    coord_t x0, coord_t y0,
    coord_t x1, coord_t y1
    ) {
    al_coord_t dx0 = static_cast<al_coord_t> (x0) - static_cast<al_coord_t> (commonX);
    al_coord_t dy0 = static_cast<al_coord_t> (y0) - static_cast<al_coord_t> (commonY);
    al_coord_t dx1 = static_cast<al_coord_t> (x1) - static_cast<al_coord_t> (commonX);
    al_coord_t dy1 = static_cast<al_coord_t> (y1) - static_cast<al_coord_t> (commonY);
    al_coord_t ret = dx0 * dy1 - dy0*dx1;
    return static_cast<coord_t> (ret);
  }

  static void products(
    coord_t commonX, coord_t commonY,
    coord_t x0, coord_t y0,
    coord_t x1, coord_t y1,
    coord_t& dotProd, coord_t& crossProd
    ) {
    al_coord_t dx0 = static_cast<al_coord_t> (x0) - static_cast<al_coord_t> (commonX);
    al_coord_t dy0 = static_cast<al_coord_t> (y0) - static_cast<al_coord_t> (commonY);
    al_coord_t dx1 = static_cast<al_coord_t> (x1) - static_cast<al_coord_t> (commonX);
    al_coord_t dy1 = static_cast<al_coord_t> (y1) - static_cast<al_coord_t> (commonY);
    al_coord_t dp = dx0 * dx1 + dy0*dy1;
    al_coord_t cp = dx0 * dy1 - dx1*dy0;
    dotProd = static_cast<coord_t> (dp);
    crossProd = static_cast<coord_t> (cp);
  }

  static void direction(
    coord_t x0, coord_t y0,
    coord_t x1, coord_t y1,
    coord_t& dx, coord_t& dy
    ) {
    dx = x1 - x0;
    dy = y1 - y0;
  }

  static void left_ortho_direction(
      coord_t x0, coord_t y0,
      coord_t x1, coord_t y1,
      coord_t& dx, coord_t&dy
  ) {
    dx = y0 - y1;
    dy = x1 - x0;
  }

  static void right_ortho_direction(
      coord_t x0, coord_t y0,
      coord_t x1, coord_t y1,
      coord_t& dx, coord_t&dy
  ) {
    dx = y1 - y0;
    dy = x0 - x1;
  }

  static bool versor(
    coord_t x0, coord_t y0,
    coord_t x1, coord_t y1,
    coord_t& dx, coord_t& dy,
    const eps_prec<al_coord_t>& prec = eps_prec<al_coord_t>::DefaultPrec
  ) {
    al_coord_t ddx = static_cast<al_coord_t> (x1) - static_cast<al_coord_t> (x0);
    al_coord_t ddy = static_cast<al_coord_t> (y1) - static_cast<al_coord_t> (y0);
    al_coord_t norm = std::hypot(ddx, ddy);
    bool ret = norm > prec.eps();
    if(ret) {
      dx=ddx/norm; dy=ddy/norm;
    }
    return ret;
  }

  // the line is defined by start point and direction (dir=endPoint-startPoint)
  // Returns true if the {x,y} point lays on the line, in which case lineP is
  //  updated with the t value that makes p=={anchorX+t*dirX,anchorY+t*dirY}
  static bool p_on_dirline_param(
    coord_t anchorX, coord_t anchorY,
    coord_t dirX, coord_t dirY,
    coord_t x, coord_t y,
    al_coord_t& lineP,
    linear_variety lineType=linear_variety::line,
    const eps_prec<al_coord_t>& prec = eps_prec<al_coord_t>::DefaultPrec
  )
  {
    al_coord_t param=0;
    bool ret=p_on_dirline(
      static_cast<al_coord_t>(anchorX), static_cast<al_coord_t>(anchorY),
      static_cast<al_coord_t>(dirX), static_cast<al_coord_t>(dirY),
      static_cast<al_coord_t>(x), static_cast<al_coord_t>(y),
      param, lineType, prec
    );
    if(ret) lineP=static_cast<coord_t>(param);
    return ret;
  }


  static bool intersection_params(
    coord_t s0x, coord_t s0y, coord_t e0x, coord_t e0y,
    coord_t s1x, coord_t s1y, coord_t e1x, coord_t e1y,
    al_coord_t& l0param, al_coord_t& l1param,
    linear_variety l0type=linear_variety::line,
    linear_variety l1type=linear_variety::line,
    const eps_prec<al_coord_t>& prec = eps_prec<al_coord_t>::DefaultPrec
  )
  {
    bool ret=false;
    al_coord_t d0x = static_cast<al_coord_t> (e0x) - static_cast<al_coord_t> (s0x);
    al_coord_t d0y = static_cast<al_coord_t> (e0y) - static_cast<al_coord_t> (s0y);
    al_coord_t d1x = static_cast<al_coord_t> (e1x) - static_cast<al_coord_t> (s1x);
    al_coord_t d1y = static_cast<al_coord_t> (e1y) - static_cast<al_coord_t> (s1y);
    bool seg0Degenerated = (d0x * d0x + d0y * d0y) <= prec.sqeps();
    bool seg1Degenerated = (d1x * d1x + d1y * d1y) <= prec.sqeps();
    if (!seg0Degenerated && !seg1Degenerated) {
      al_coord_t det = d0x * d1y - d1x*d0y;
      if (std::abs(det) <= prec.sqeps()) { // parallel or coincident
        al_coord_t p1s_param=0, p1e_param=0;
        ret=p_on_dirline(
          static_cast<al_coord_t> (s0x), static_cast<al_coord_t> (s0y),
          d0x, d0y,
          static_cast<al_coord_t> (s1x), static_cast<al_coord_t> (s1y),
          p1s_param, linear_variety::line, prec
        );
        if(ret) {
          // See how the end of second line/ray/segment is positioned
          ret=p_on_dirline(
            static_cast<al_coord_t> (s0x), static_cast<al_coord_t> (s0y),
            d0x, d0y,
            static_cast<al_coord_t> (e1x), static_cast<al_coord_t> (e1y),
            p1e_param, linear_variety::line, prec
          );
        }
        if(ret) { // colinear segs/rays/lines
          switch(l0type) {
            case linear_variety::line:
              // ret is true already, we'll consider the "mid point" of the second line
              // as the intersection.
              l0param=(p1s_param+p1e_param)/2;
              l1param=0.5;
              break;
            case linear_variety::ray:
              switch(l1type) {
                case linear_variety::line:
                  // the second line intersects with the first ray anyway
                  if(p1s_param<0 && p1e_param<0) {
                    // second line is defined by points behind the start of first ray
                    ret=l0param=0;
                    p_on_dirline(
                      static_cast<al_coord_t> (s1x), static_cast<al_coord_t> (s1y),
                      d1x, d1y,
                      static_cast<al_coord_t> (s0x), static_cast<al_coord_t> (s0y),
                      l1param, linear_variety::line, prec
                    );
                  }
                  else if(p1s_param>=0 && p1e_param>=0) {
                    // take the "mid-point" of the second line as the intersection
                    l0param=(p1s_param+p1e_param)/2;
                    l1param=0.5;
                  }
                  else if(p1s_param>=0) { // this means p1e_param <0
                    l0param=p1s_param;
                    l1param=0.0;
                  }
                  else { // p1e_param >= 0  and p1s_param <0
                    l0param=p1e_param;
                    l1param=1.0;
                  }
                  break;
                case linear_variety::ray: // collinear ray-ray intersection
                  if(p1s_param>=0 && p1e_param>=0) {
                    // first ray catches up the second one
                    // take the "mid-point" of the second ray as the intersection
                    l0param=(p1s_param+p1e_param)/2;
                    l1param=0.5;
                  }
                  else if(p1s_param>=0) { // opposite rays with intersection
                    l0param=p1s_param;
                    l1param=0.0;
                  }
                  else if(p1e_param>=0) {
                    // same dir rays, with the second end point after the first ray's source
                    l0param=p1e_param;
                    l1param=1.0;
                  }
                  else if(p1s_param<p1e_param) {
                    // the second ray is defined by points behind the first ray's source
                    // but the second ray travels towards the first one
                    l0param=0.0;
                    p_on_dirline(
                      static_cast<al_coord_t> (s1x), static_cast<al_coord_t> (s1y),
                      d1x, d1y,
                      (static_cast<al_coord_t> (s0x)+static_cast<al_coord_t> (e0x))/2,
                      (static_cast<al_coord_t> (s0y)+static_cast<al_coord_t> (e0y))/2,
                      l1param, linear_variety::ray, prec
                    );
                  }
                  else { // opposite rays not sharing any points
                    ret=false;
                  }
                  break;
                case linear_variety::segment: // collinear ray-segment intersection
                  if(p1s_param>=0 && p1e_param>=0) {
                    // second segment fully on the first ray
                    // take the "mid-point" of the second segment as the intersection
                    l0param=(p1s_param+p1e_param)/2;
                    l1param=0.5;
                  }
                  else if(p1s_param>=0) { // only the start of second segment is on the ray
                    l0param=p1s_param;
                    l1param=0.0;
                  }
                  else if(p1e_param>=0) {
                    // only the end of the second segment is on the ray
                    l0param=p1e_param;
                    l1param=1.0;
                  }
                  else { // no ends of the segment are on the actual ray
                    ret=false;
                  }
                  break;
              }
              break;
            case linear_variety::segment:
            default:
              switch(l1type) {
                case linear_variety::line:
                  // take the mid-point of first segment as the intersection
                  l0param=0.5;
                  p_on_dirline(
                    static_cast<al_coord_t> (s1x), static_cast<al_coord_t> (s1y),
                    d1x, d1y,
                    (static_cast<al_coord_t> (s0x)+static_cast<al_coord_t> (e0x))/2,
                    (static_cast<al_coord_t> (s0y)+static_cast<al_coord_t> (e0y))/2,
                    l1param, linear_variety::line, prec
                  );
                  break;
                case linear_variety::ray: // segment-ray intersection
                  if(
                       (p1s_param<=0 && p1e_param>p1s_param) // ray same direction with segm
                    || (p1s_param>=1 && p1e_param<p1s_param) // ray opposite segm
                  ) {
                    // the ray start is outside the segment, but the ray travels towards segment
                    // the entire first segment lays on the ray, take its mid as intersection
                    ret=p_on_dirline(
                      static_cast<al_coord_t> (s1x), static_cast<al_coord_t> (s1y),
                      d1x, d1y,
                      (static_cast<al_coord_t> (s0x)+static_cast<al_coord_t> (e0x))/2,
                      (static_cast<al_coord_t> (s0y)+static_cast<al_coord_t> (e0y))/2,
                      l1param, linear_variety::ray, prec
                    );
                    l0param=0.5;
                  }
                  else if(
                       (p1s_param<0 && p1e_param<p1s_param) // ray opposite segm
                    || (p1s_param>1 && p1e_param>p1s_param) // ray same direction with segm
                  ) { // the ray start is outside the segment and the ray lets the segment behind
                    // no segm on ray
                    ret=false;
                  }
                  // it can only be that the 2nd ray starts inside the 1st segment
                  else if(p1s_param<p1e_param) { // ray travels towards the segment's end
                    l0param=1.0;
                    p_on_dirline(
                      static_cast<al_coord_t> (s1x), static_cast<al_coord_t> (s1y),
                      d1x, d1y,
                      static_cast<al_coord_t> (e0x),  static_cast<al_coord_t> (e0y),
                      l1param, linear_variety::ray, prec
                    );
                  }
                  else if(p1s_param>p1e_param){ // ray travels towards the segment start
                    l0param=0.0;
                    ret=p_on_dirline(
                      static_cast<al_coord_t> (s1x), static_cast<al_coord_t> (s1y),
                      d1x, d1y,
                      static_cast<al_coord_t> (s0x),  static_cast<al_coord_t> (s0y),
                      l1param, linear_variety::ray, prec
                    );
                  }
                  else {
                    assert("Should not see this: collinear segment-ray intersection");
                  }
                  break;
                case linear_variety::segment: // segment-segment intersection
                default:  {
                  double p1min=std::min(p1s_param, p1e_param);
                  double p1max=std::max(p1s_param, p1e_param);
                  if(p1min>=0 && p1max<=1) {
                    // second segment embedded in the first. Take it's mid point
                    l0param=(p1s_param+p1e_param)/2.0;
                    l1param=0.5;
                  }
                  else if(p1min<=0 && p1max>=1) {
                    // first segment embedded in the second, take it's mid point
                    l0param=0.5;
                    ret=p_on_dirline(
                      static_cast<al_coord_t> (s1x), static_cast<al_coord_t> (s1y),
                      d1x, d1y,
                      (static_cast<al_coord_t> (s0x)+static_cast<al_coord_t> (e0x))/2.0,
                      (static_cast<al_coord_t> (s0y)+static_cast<al_coord_t> (e0y))/2.0,
                      l1param, linear_variety::ray, prec
                    );
                  }
                  else if(
                       (p1s_param<0 && p1e_param<0)
                    || (p1s_param>1.0 && p1e_param>1.0)
                  ) { // no overlapping segments
                    ret=false;
                  }
                  // one end of the second segment is inside, the other is outside the first segment
                  else if(p1s_param>=0 && p1s_param<=1) {
                    // second segment has its start inside the first. The end of the second must be outside
                    l0param=
                        (p1s_param<p1e_param)
                      ? l0param=(p1s_param+1.0)/2.0 // the end of the first segment is inside the second segment
                      : p1s_param/2.0 // the start of the first segment is inside the second
                    ;
                    ret=p_on_dirline(
                      static_cast<al_coord_t> (s1x), static_cast<al_coord_t> (s1y),
                      d1x, d1y,
                      static_cast<al_coord_t> (s0x)+d0x*l0param,
                      static_cast<al_coord_t> (s0y)+d0y*l0param,
                      l1param, linear_variety::segment, prec
                    );
                  }
                  else if(p1e_param>=0 && p1e_param<=1) {
                    // the end of the second segment is inside the first
                    l0param=
                        (p1e_param<p1s_param)
                      ? (p1e_param+1.0)/2.0
                      : p1e_param/2.0;
                    ret=p_on_dirline(
                      static_cast<al_coord_t> (s1x), static_cast<al_coord_t> (s1y),
                      d1x, d1y,
                      static_cast<al_coord_t> (s0x)+d0x*l0param,
                      static_cast<al_coord_t> (s0y)+d0y*l0param,
                      l1param, linear_variety::segment, prec
                    );
                  }
                  else {
                    assert("Should not see this: collinear segment segment intersection");
                  }
                }
                break;
              }
              break;
          }
        }
      }
      else { // non-degenerated with a significant determinant
        al_coord_t dsx=static_cast<al_coord_t>(s1x)-static_cast<al_coord_t>(s0x);
        al_coord_t dsy=static_cast<al_coord_t>(s1y)-static_cast<al_coord_t>(s0y);
        al_coord_t param0=(dsx*d1y-dsy*d1x)/det;
        al_coord_t param1=(dsx*d0y-dsy*d0x)/det;
        switch(l0type) {
          case linear_variety::segment:
            ret = (param0>=0 && param0<=1);
            break;
          case linear_variety::ray:
            ret = (param0>=0);
            break;
          case linear_variety::line:
          default:
            ret=true;
            break;
        }
        if(ret) {
          switch(l1type) {
            case linear_variety::segment:
              ret=(param1>=0 && param1<=1);
              break;
            case linear_variety::ray:
              ret=(param1>=0);
              break;
            case linear_variety::line:
            default:
              ret=true;
              break;
          }
        }
        if(ret) {
          l0param=param0;
          l1param=param1;
        }
      }
    }
    else if(seg0Degenerated && seg0Degenerated) {
      // two degenerated lines (two points actually) "intersect"
      // only if their "midpoints" coincide
      al_coord_t m0x=(static_cast<al_coord_t>(s0x)+static_cast<al_coord_t>(e0x))/2;
      al_coord_t m0y=(static_cast<al_coord_t>(s0y)+static_cast<al_coord_t>(e0y))/2;
      al_coord_t m1x=(static_cast<al_coord_t>(s1x)+static_cast<al_coord_t>(e1x))/2;
      al_coord_t m1y=(static_cast<al_coord_t>(s1y)+static_cast<al_coord_t>(e1y))/2;
      m1x-=m0x;
      m1y-=m0y;
      ret=(m1x*m1x+m1y*m1y) < prec.sqeps();
      if(ret) {
        l0param=l1param=0.5;
      }
    }
    else if(seg0Degenerated) {
      // will "intersect" of if the midpoint of seg 0 lays onto the line1
      al_coord_t m0x=(static_cast<al_coord_t>(s0x)+static_cast<al_coord_t>(e0x))/2;
      al_coord_t m0y=(static_cast<al_coord_t>(s0y)+static_cast<al_coord_t>(e0y))/2;
      al_coord_t p=0;
      ret=p_on_dirline(
        static_cast<al_coord_t>(s1x), static_cast<al_coord_t>(s1y),
        d1x, d1y,
        m0x, m0y,
        p, l1type, prec
      );
      if(ret) {
        l0param=0.5;
        l1param=p;
      }
    }
    else { // seg1 degenerated
      // will "intersect" of if the midpoint of seg 0 lays onto the line1
      al_coord_t m1x=(static_cast<al_coord_t>(s1x)+static_cast<al_coord_t>(e1x))/2;
      al_coord_t m1y=(static_cast<al_coord_t>(s1y)+static_cast<al_coord_t>(e1y))/2;
      al_coord_t p=0;
      ret=p_on_dirline(
        static_cast<al_coord_t>(s0x), static_cast<al_coord_t>(s0y),
        d0x, d0y,
        m1x, m1y,
        p, l1type, prec
      );
      if(ret) {
        l1param=0.5;
        l0param=p;
      }
    }
    return ret;
  }
};

} // namespace geomalgos2d
#endif /* LINEAR_GEOM_HPP */
