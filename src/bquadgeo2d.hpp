#ifndef BQUADGEO2D_HPP
#define BQUADGEO2D_HPP

#include "linear2d.hpp"

namespace geomalgos2d {

template <typename coord_t, typename al_coord_t>
coord_t quadcomb(coord_t start, coord_t ctrl, coord_t end, al_coord_t param)
{
  al_coord_t s = static_cast<al_coord_t>(start);
  al_coord_t c = static_cast<al_coord_t>(ctrl);
  al_coord_t e = static_cast<al_coord_t>(end);
  al_coord_t t = static_cast<al_coord_t>(param);
  al_coord_t onet = static_cast<al_coord_t>(1) - t;

  al_coord_t ret = onet*(onet*s + t*c) + t*(onet*c + t*e);
  return static_cast<coord_t>(ret);
}

template <typename coord_t>
coord_t quadcomb(coord_t start, coord_t ctrl, coord_t end, coord_t param)
{
  coord_t onet = param;

  coord_t ret = onet(onet*start + param*ctrl) + param*(onet*ctrl + param*end);
  return ret;
}


template <typename coord_t, typename al_coord_t = double> class bquadgeo2d {
  static bool quadvals_internal(
    al_coord_t sx, al_coord_t sy,
    al_coord_t cx, al_coord_t cy,
    al_coord_t ex, al_coord_t ey,
    al_coord_t param,
    al_coord_t* px = nullptr, al_coord_t* py = nullptr,
    al_coord_t* vdirx = nullptr, al_coord_t* vdiry = nullptr,
    al_coord_t* adirx = nullptr, al_coord_t* adiry = nullptr
  ) {
    bool needpos = px || py;
    bool needv = vdirx || vdiry;
    bool needa = adirx || adiry;
    if(needpos || needv || needa) {
      al_coord_t onet = static_cast<al_coord_t>(1) - param;
      if(px) *px = onet*(onet*sx + param*cx) + param*(onet*cx + param*ex);
      if(py) *py = onet*(onet*sy + param*cy) + param*(onet*cy + param*ey);
      if(needv || needa) {
        al_coord_t d0x = cx - sx, d0y = cy - sy;
        al_coord_t d1x = ex - cx, d1y = ey - cy;
        if(vdirx) *vdirx = 2*(onet*d0x + param*d1x);
        if(vdiry) *vdiry = 2*(onet*d0y + param*d1y);
        if(adirx) *adirx = 2*(d1x - d0x);
        if(adiry) *adiry = 2*(d1y - d0y);
      }
    }
    return true;
  }

public:
  static bool quadvals(
    coord_t startx, coord_t starty,
    coord_t ctrlx, coord_t ctrly,
    coord_t endx, coord_t endy,
    coord_t param,
    coord_t* posx = nullptr, coord_t* posy = nullptr,
    coord_t* vdirx = nullptr, coord_t* vdiry = nullptr,
    coord_t* adirx = nullptr, coord_t* adiry = nullptr
  ) {
    bool needpos = posx || posy;
    bool needv = vdirx || vdiry;
    bool needa = adirx || adiry;
    if(needpos || needv || needa) {
      al_coord_t px, py, vx, vy, ax, ay;
      quadvals_internal(
        static_cast<al_coord_t>(startx), static_cast<al_coord_t>(starty),
        static_cast<al_coord_t>(ctrlx), static_cast<al_coord_t>(ctrly),
        static_cast<al_coord_t>(endx), static_cast<al_coord_t>(endy),
        static_cast<al_coord_t>(param),
        posx ? &px : nullptr, posy ? &py : nullptr,
        vdirx ? &vx : nullptr, vdiry ? &vy : nullptr,
        adirx ? &ax : nullptr, adiry ? &ay : nullptr
      );
      if(posx) *posx = static_cast<coord_t>(px);
      if(posy) *posy = static_cast<coord_t>(py);
      if(vdirx) *vdirx = static_cast<coord_t>(vx);
      if(vdiry) *vdiry = static_cast<coord_t>(vy);
      if(adirx) *adirx = static_cast<coord_t>(ax);
      if(adiry) *adiry = static_cast<coord_t>(ay);
    }
    return true;
  }

  static bool quadprops(
    coord_t sx, coord_t sy,
    coord_t cx, coord_t cy,
    coord_t ex, coord_t ey,
    coord_t* vertexParam = nullptr,
    coord_t* vertexx = nullptr, coord_t* vertexy = nullptr,
    coord_t* focusx = nullptr, coord_t* focusy = nullptr,
    const eps_prec<al_coord_t>& prec = eps_prec<al_coord_t>::DefaultPrec
  ) {
    bool ret = false;
    if(vertexParam || vertexx || vertexy || focusx || focusy) {
      al_coord_t vx, vy, ax, ay;
      quadvals_internal(
        sx, sy, cx, cy, ex, ey,
        static_cast<al_coord_t>(0),
        nullptr, nullptr,
        &vx, &vy, &ax, &ay
      );
      al_coord_t amodulus = ax*ax + ay*ay;
      ret = amodulus > prec.sqeps();
      if(ret) { // not a line travelled at constant speed
        al_coord_t px, py; // vertex coords here, if needed
        al_coord_t vt = -dotp<al_coord_t>(vx, vy, ax, ay)/amodulus;
        if(vertexParam) *vertexParam = static_cast<coord_t>(vt);
        bool needFocus = focusx || focusy;
        if(vertexx || vertexy || needFocus) {
          quadvals_internal(
            sx, sy, cx, cy, ex, ey,
            vt,
            vertexx || needFocus ? &px : nullptr, vertexy || needFocus ? &py : nullptr,
            needFocus ? &vx : nullptr, needFocus ? &vy : nullptr // acceleration stays the same, don't ask for it again
          );
          if(vertexx) *vertexx = static_cast<coord_t>(px);
          if(vertexy) *vertexy = static_cast<coord_t>(py);
          if(needFocus) {
            // Focus is midpoint betweem vertex and the center of the osculating circle tangent in the vertex
            // see https://mathworld.wolfram.com/OsculatingCircle.html

            al_coord_t vcrossa = crossp(vx, vy, ax, ay);
            ret = std::abs(vcrossa) > prec.sqeps();
            if(ret) { // quad not degenerated in a line travelled at constant acceleration
              al_coord_t coef=(vx*vx+vy*vy)/vcrossa;
              if(focusx) *focusx = static_cast<coord_t>(px - coef*vy/static_cast<al_coord_t>(2));
              if(focusy) *focusy = static_cast<coord_t>(py + coef*vx/static_cast<al_coord_t>(2));
            }
            // else degenerated in a line
          } // else don't need focus position
        } // else don't need vertex nor focus
      } // else the parabola is totally degenerated to a line travelled at constant speed
    } // else nothing was requested to compute
    return ret;
  }
};

} // namespace geomalgos2d {

#endif // BQUADGEO2D_HPP
