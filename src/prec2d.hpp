#ifndef PREC2D_HPP
#define PREC2D_HPP

#include <limits>
#include <cmath>

#include "dist2d.hpp"

namespace geomalgos2d
{
template <typename al_coord_t> class eps_prec {
private:
  al_coord_t m_SqEps;
  al_coord_t m_Eps;
public:
  static const eps_prec DefaultPrec;

  eps_prec() { eps(-1); }
  eps_prec(al_coord_t v) {eps(v);}

  al_coord_t eps() const { return m_Eps; }

  al_coord_t sqeps() const { return m_SqEps; }

  al_coord_t eps(al_coord_t newVal=-1) {
    al_coord_t ret=m_Eps;
    if(newVal<0) {
      al_coord_t e=std::numeric_limits<al_coord_t>::epsilon()*1024.0;
      newVal=std::sqrt(e);
    }
    m_Eps=newVal;
    m_SqEps=newVal*newVal;
    return ret;
  }

  template <
      typename coord_t = al_coord_t,
      coord_t (*dist_func)(coord_t, coord_t, coord_t, coord_t) = dist2d<coord_t, al_coord_t>::euclid
  >
  static bool same_point_test(
    coord_t x0, coord_t y0, coord_t x1, coord_t y1,
    const eps_prec<al_coord_t>& prec = DefaultPrec
  ) {
    return dist_func(x0, y0, x1, y1) <= prec.eps();
  }

  template <
      typename coord_t,
      coord_t (*dist_func_sq)(coord_t, coord_t, coord_t, coord_t) = dist2d<coord_t, al_coord_t>::sqeuclid
  >
  static bool same_point_test_sq(
    coord_t x0, coord_t y0, coord_t x1, coord_t y1,
    const eps_prec<al_coord_t>& prec = DefaultPrec
  ) {
    return dist_func_sq(x0, y0, x1, y1) <= prec.sqeps();
  }

};

template <typename al_coord_t> const eps_prec<al_coord_t> eps_prec<al_coord_t>::DefaultPrec =
    eps_prec<al_coord_t>(-1);

} // namespace geomalgos2d

#endif // PREC2D_HPP
