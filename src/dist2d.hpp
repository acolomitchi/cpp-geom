#ifndef DIST2D_HPP
#define DIST2D_HPP

#include <cmath>

namespace geomalgos2d
{

template <typename coord_t, typename al_coord_t> class dist2d {
public:
  static coord_t sqeuclid(coord_t x0, coord_t y0, coord_t x1, coord_t y1){
    al_coord_t dx = static_cast<al_coord_t> (x1) - static_cast<al_coord_t> (x0);
    al_coord_t dy = static_cast<al_coord_t> (y1) - static_cast<al_coord_t> (y0);
    return static_cast<coord_t> (dx * dx + dy * dy);
  }

  static coord_t euclid(coord_t x0, coord_t y0, coord_t x1, coord_t y1) {
    al_coord_t dx = static_cast<al_coord_t> (x1) - static_cast<al_coord_t> (x0);
    al_coord_t dy = static_cast<al_coord_t> (y1) - static_cast<al_coord_t> (y0);
    return static_cast<coord_t> (std::hypot(dx, dy));
  }
};
}
#endif // DIST2D_HPP
