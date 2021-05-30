#include <catch2/catch.hpp>

#include <limits>

#include "sorta_quad.hpp"
#define NUM_PI 3.141592653589793238462643383279502884

using quad = geomalgos2d::bquadgeo2d<double, double>;


TEST_CASE("basic quad vals", "[quad]") {
  sorta_quad qb(-1, 1, 0, -1, 1, 1);
  SECTION("vals trivial") {
    sorta_quad q(qb);
    double px, py, vx, vy, ax, ay;
    bool r = quad::quadvals(q.sx, q.sy, q.cx, q.cy, q.ex, q.ey, 0.5, &px, &py, &vx, &vy, &ax, &ay);
    REQUIRE(r);
    REQUIRE(px == 0); REQUIRE(py == 0);
    REQUIRE(vx == 2); REQUIRE(vy == 0);
    REQUIRE(ax == 0); REQUIRE(ay == 8);
  }
  SECTION("vals non-trivial") {
    using prec = geomalgos2d::eps_prec<double>;
    sorta_quad q = qb.rotated(0.7524986583);
    double px, py, vx, vy, ax, ay;
    bool r = quad::quadvals(q.sx, q.sy, q.cx, q.cy, q.ex, q.ey, 0.5, &px, &py, &vx, &vy, &ax, &ay);
    REQUIRE(r);
    REQUIRE(px == 0); REQUIRE(py == 0);
    double v = std::hypot(vx, vy);
    REQUIRE(std::abs(v-2.0) < prec::DefaultPrec.eps());
    double a = std::hypot(ax, ay);
    REQUIRE(std::abs(a - 8.0) < prec::DefaultPrec.eps());
  }
  SECTION("Null ptrs on y") {
    sorta_quad q = qb.rotated(0.7524986583);
    double px, vx, ax;
    bool r = quad::quadvals(q.sx, q.sy, q.cx, q.cy, q.ex, q.ey, 0.5, &px, nullptr, &vx, nullptr, &ax, nullptr);
    REQUIRE(r);
  }
  SECTION("Null ptrs on x") {
    sorta_quad q = qb.rotated(0.7524986583);
    double py, vy, ay;
    bool r = quad::quadvals(q.sx, q.sy, q.cx, q.cy, q.ex, q.ey, 0.5, nullptr, &py, nullptr, &vy, nullptr, &ay);
    REQUIRE(r);
  }
}

TEST_CASE("basic quad props", "[quad]") {
  sorta_quad qb(-1, 1, 0, 0, 1, 1);
  SECTION("props trivial") {
    sorta_quad q(qb);
    double t, vx, vy, fx, fy;
    bool r = quad::quadprops(
          q.sx, q.sy, q.cx, q.cy, q.ex, q.ey,
          &t, &vx, &vy, &fx, &fy);
    REQUIRE(r);
    REQUIRE(t == 0.5);
    REQUIRE(vx == 0); REQUIRE(vy == 0.5);
    REQUIRE(fx == 0); REQUIRE(fy == 1);
  }
  SECTION("props non-trivial") {
    using prec = geomalgos2d::eps_prec<double>;
    sorta_quad q = qb.translated(0, -1);
    q.rotate(0.86415469792645);
    q.translate(0, 1);
    double t, vx, vy, fx, fy;
    bool r = quad::quadprops(
          q.sx, q.sy, q.cx, q.cy, q.ex, q.ey,
          &t, &vx, &vy, &fx, &fy);
    REQUIRE(r);
    REQUIRE(t == 0.5);
    // REQUIRE(vx == 0); REQUIRE(vy == 0.5);
    REQUIRE(std::abs(fx) < prec::DefaultPrec.eps());
    REQUIRE(std::abs(fy - 1) < prec::DefaultPrec.eps());
  }
}

