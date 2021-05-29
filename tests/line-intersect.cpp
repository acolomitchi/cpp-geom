#include <catch2/catch.hpp>

#include <limits>
#include "sorta_line.hpp"

using namespace geomalgos2d;

TEST_CASE("line intersection trivial", "[positive][line]") {
  sorta_line l0(0, 0, 1, 1, linear_variety::line);
  SECTION("with line") {
    sorta_line l1(0, 1, 1, 0, linear_variety::line);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("with ray") {
    sorta_line l1(0, 1, 1, 0, linear_variety::ray);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("with segment") {
    sorta_line l1(0, 1, 1, 0, linear_variety::segment);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
}

TEST_CASE("line intersection same support", "[positive][line]") {
  sorta_line l0(0, 0, 1, 1, linear_variety::line);
  SECTION("with line") {
    sorta_line l1(0.25, 0.25, 1, 1, linear_variety::line);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("with ray") {
    sorta_line l1(0.25, 0.25, 1, 1, linear_variety::ray);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("with segment") {
    sorta_line l1(0.25, 0.25, 1, 1, linear_variety::segment);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
}

TEST_CASE("line intersection through end ray/segments", "[positive][line]") {
  sorta_line l0(0, 0, 1, 1, linear_variety::line);
  SECTION("with ray start") {
    sorta_line l1(0.5, 0.5, 1, 1, linear_variety::ray);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("with segment begin") {
    sorta_line l1(0.5, 0.5, 1, 1, linear_variety::segment);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("with segment end") {
    sorta_line l1(1, 1, 0.5, 0.5, linear_variety::segment);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
}

TEST_CASE("line intersection negative parallel", "[negative][line][ray][segment]") {
  sorta_line l0(0, 0, 1, 1, linear_variety::line);
  SECTION("parallel lines") {
    sorta_line l1(0, -1, 1, 0, linear_variety::line);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE_FALSE(intersectionExists);
    REQUIRE(l0p == std::numeric_limits<double>::max());
    REQUIRE(l1p == std::numeric_limits<double>::max());
    REQUIRE(l0x == std::numeric_limits<double>::min());
    REQUIRE(l0y == std::numeric_limits<double>::min());
    REQUIRE(l1x == std::numeric_limits<double>::min());
    REQUIRE(l1y == std::numeric_limits<double>::min());
  }
  SECTION("anti-parallel lines") {
    sorta_line l1(1, 0, 0, -1, linear_variety::line);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE_FALSE(intersectionExists);
    REQUIRE(l0p == std::numeric_limits<double>::max());
    REQUIRE(l1p == std::numeric_limits<double>::max());
    REQUIRE(l0x == std::numeric_limits<double>::min());
    REQUIRE(l0y == std::numeric_limits<double>::min());
    REQUIRE(l1x == std::numeric_limits<double>::min());
    REQUIRE(l1y == std::numeric_limits<double>::min());
  }
  SECTION("with parallel ray") {
    sorta_line l1(0, -1, 1, 0, linear_variety::ray);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE_FALSE(intersectionExists);
    REQUIRE(l0p == std::numeric_limits<double>::max());
    REQUIRE(l1p == std::numeric_limits<double>::max());
    REQUIRE(l0x == std::numeric_limits<double>::min());
    REQUIRE(l0y == std::numeric_limits<double>::min());
    REQUIRE(l1x == std::numeric_limits<double>::min());
    REQUIRE(l1y == std::numeric_limits<double>::min());
  }
  SECTION("with anti-parallel ray") {
    sorta_line l1(1, 0, 0, -1, linear_variety::ray);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE_FALSE(intersectionExists);
    REQUIRE(l0p == std::numeric_limits<double>::max());
    REQUIRE(l1p == std::numeric_limits<double>::max());
    REQUIRE(l0x == std::numeric_limits<double>::min());
    REQUIRE(l0y == std::numeric_limits<double>::min());
    REQUIRE(l1x == std::numeric_limits<double>::min());
    REQUIRE(l1y == std::numeric_limits<double>::min());
  }
  SECTION("with parallel segment") {
    sorta_line l1(0, -1, 1, 0, linear_variety::segment);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE_FALSE(intersectionExists);
    REQUIRE(l0p == std::numeric_limits<double>::max());
    REQUIRE(l1p == std::numeric_limits<double>::max());
    REQUIRE(l0x == std::numeric_limits<double>::min());
    REQUIRE(l0y == std::numeric_limits<double>::min());
    REQUIRE(l1x == std::numeric_limits<double>::min());
    REQUIRE(l1y == std::numeric_limits<double>::min());
  }
  SECTION("with anti-parallel segment") {
    sorta_line l1(1, 0, 0, -1, linear_variety::segment);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE_FALSE(intersectionExists);
    REQUIRE(l0p == std::numeric_limits<double>::max());
    REQUIRE(l1p == std::numeric_limits<double>::max());
    REQUIRE(l0x == std::numeric_limits<double>::min());
    REQUIRE(l0y == std::numeric_limits<double>::min());
    REQUIRE(l1x == std::numeric_limits<double>::min());
    REQUIRE(l1y == std::numeric_limits<double>::min());
  }
}

TEST_CASE("ray/(ray,segment) intersection in one end", "[positive][ray][segment") {
  sorta_line l0(0, -1, 1, 0, linear_variety::ray);
  SECTION("rays sharing start, different directions") {
    sorta_line l1(0, -1, 0, 1, linear_variety::ray);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(l0p == 0);
    REQUIRE(l1p == 0);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("rays sharing start, opposite directions") {
    sorta_line l1(0, -1, -1, 2, linear_variety::ray);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(l0p == 0);
    REQUIRE(l1p == 0);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("ray/segment sharing start, different directions") {
    sorta_line l1(0, -1, 0, 1, linear_variety::segment);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(l0p == 0);
    REQUIRE(l1p == 0);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("ray/segment sharing end, different directions") {
    sorta_line l1(0, 1, 0, -1, linear_variety::segment);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(l0p == 0);
    REQUIRE(l1p == 1);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("ray/segment sharing start, opposite directions") {
    sorta_line l1(0, -1, -1, 2, linear_variety::segment);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(l0p == 0);
    REQUIRE(l1p == 0);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("ray/segment sharing segment end, opposite directions") {
    sorta_line l1(-1, 2, 0, -1, linear_variety::segment);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(l0p == 0);
    REQUIRE(l1p == 1);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("ray/segment with segment start on ray") {
    sorta_line l1(0.5, -0.5, -1, 2, linear_variety::segment);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(l1p == 0);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("ray/segment with segment end on ray") {
    sorta_line l1(-1, 2, 0.5, -0.5, linear_variety::segment);
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(l1p == 1);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
}

TEST_CASE("same dir ray/(ray,segment) intersection", "[positive][ray][segment") {
  sorta_line l0(0, -1, 1, 0, linear_variety::ray);
  SECTION("with ovelaping ray downward of original") {
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    sorta_line l1(2, 1, 3, 2, linear_variety::ray);
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("with ovelaping ray upward of original") {
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    sorta_line l1(-1, -2, 0, -1, linear_variety::ray);
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("with opposite ray") {
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    sorta_line l1(2, 1, 0.5, -0.5, linear_variety::ray);
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("with same dir ovelaping segment downward of original") {
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    sorta_line l1(-2,-3, 0.5, -0.5, linear_variety::segment);
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("with opposite dir ovelaping segment downward of original") {
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    sorta_line l1(0.5, -0.5, -2,-3, linear_variety::segment);
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("with segment totally embedded in the ray") {
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    sorta_line l1(3, 2, 5, 4, linear_variety::ray);
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
}

TEST_CASE("ray/ray negative parallel", "[negative][ray]") {
  sorta_line l0(0, -1, 1, 0, linear_variety::ray);
  SECTION("parallel ray") {
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    sorta_line l1(0, 1, 1, 2, linear_variety::ray);
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE_FALSE(intersectionExists);
    REQUIRE(l0p == std::numeric_limits<double>::max());
    REQUIRE(l1p == std::numeric_limits<double>::max());
    REQUIRE(l0x == std::numeric_limits<double>::min());
    REQUIRE(l0y == std::numeric_limits<double>::min());
    REQUIRE(l1x == std::numeric_limits<double>::min());
    REQUIRE(l1y == std::numeric_limits<double>::min());
  }
  SECTION("anti-parallel ray") {
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    sorta_line l1(1, 2, 0, 1, linear_variety::ray);
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE_FALSE(intersectionExists);
    REQUIRE(l0p == std::numeric_limits<double>::max());
    REQUIRE(l1p == std::numeric_limits<double>::max());
    REQUIRE(l0x == std::numeric_limits<double>::min());
    REQUIRE(l0y == std::numeric_limits<double>::min());
    REQUIRE(l1x == std::numeric_limits<double>::min());
    REQUIRE(l1y == std::numeric_limits<double>::min());
  }
  SECTION("parallel segment") {
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    sorta_line l1(0, 1, 1, 2, linear_variety::segment);
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE_FALSE(intersectionExists);
    REQUIRE(l0p == std::numeric_limits<double>::max());
    REQUIRE(l1p == std::numeric_limits<double>::max());
    REQUIRE(l0x == std::numeric_limits<double>::min());
    REQUIRE(l0y == std::numeric_limits<double>::min());
    REQUIRE(l1x == std::numeric_limits<double>::min());
    REQUIRE(l1y == std::numeric_limits<double>::min());
  }
  SECTION("anti-parallel segment") {
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    sorta_line l1(1, 2, 0, 1, linear_variety::segment);
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE_FALSE(intersectionExists);
    REQUIRE(l0p == std::numeric_limits<double>::max());
    REQUIRE(l1p == std::numeric_limits<double>::max());
    REQUIRE(l0x == std::numeric_limits<double>::min());
    REQUIRE(l0y == std::numeric_limits<double>::min());
    REQUIRE(l1x == std::numeric_limits<double>::min());
    REQUIRE(l1y == std::numeric_limits<double>::min());
  }
}

TEST_CASE("ray/segment with no intersection", "[negative][ray][segment") {
  sorta_line l0(0, -1, 1, 0, linear_variety::ray);
  SECTION("segment with same support downward of original") {
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    sorta_line l1(-2, -3, -1, -2, linear_variety::segment);
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE_FALSE(intersectionExists);
    REQUIRE(l0p == std::numeric_limits<double>::max());
    REQUIRE(l1p == std::numeric_limits<double>::max());
    REQUIRE(l0x == std::numeric_limits<double>::min());
    REQUIRE(l0y == std::numeric_limits<double>::min());
    REQUIRE(l1x == std::numeric_limits<double>::min());
    REQUIRE(l1y == std::numeric_limits<double>::min());
  }
  SECTION("segment with intersecting direction downwards of ray") {
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    sorta_line l1(-2, -0.5, 2, -3, linear_variety::segment);
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE_FALSE(intersectionExists);
    REQUIRE(l0p == std::numeric_limits<double>::max());
    REQUIRE(l1p == std::numeric_limits<double>::max());
    REQUIRE(l0x == std::numeric_limits<double>::min());
    REQUIRE(l0y == std::numeric_limits<double>::min());
    REQUIRE(l1x == std::numeric_limits<double>::min());
    REQUIRE(l1y == std::numeric_limits<double>::min());
  }
  SECTION("segment with intersecting direction upwards of ray origin but no intersection") {
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    sorta_line l1(-1, 1, 0, -0.5, linear_variety::segment);
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE_FALSE(intersectionExists);
    REQUIRE(l0p == std::numeric_limits<double>::max());
    REQUIRE(l1p == std::numeric_limits<double>::max());
    REQUIRE(l0x == std::numeric_limits<double>::min());
    REQUIRE(l0y == std::numeric_limits<double>::min());
    REQUIRE(l1x == std::numeric_limits<double>::min());
    REQUIRE(l1y == std::numeric_limits<double>::min());
  }
}

TEST_CASE("segment/segment intersection","[positive][segment]") {
  sorta_line l0(0, -1, 1, 0, linear_variety::segment);
  SECTION("pure inner intersection") {
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    sorta_line l1(0, 0, 1, -1, linear_variety::segment);
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("overlapping along the same direction") {
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    sorta_line l1(-1, -2, 0.5, -0.5, linear_variety::segment);
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("segment totally included in the original segment") {
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    sorta_line l1(0.25, -0.75, 0.75, -0.25, linear_variety::segment);
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("segment totally included in the original segment, shares start") {
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    sorta_line l1(0, -1, 0.75, -0.25, linear_variety::segment);
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("original segment totally included in the segment") {
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    sorta_line l1(-0.25, -1.25, 1.25, 0.25, linear_variety::segment);
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }
  SECTION("original segment totally included in the segment, shares start") {
    double l0p, l1p, l0x, l0y, l1x, l1y;
    l0p = l1p = std::numeric_limits<double>::max();
    l0x = l0y = l1x = l1y = std::numeric_limits<double>::min();
    sorta_line l1(0, -1, 1.25, 0.25, linear_variety::segment);
    bool intersectionExists=l0.intersection(l1, l0p, l1p, &l0x, &l0y, &l1x, &l1y);
    REQUIRE(intersectionExists);
    REQUIRE(eps_prec<double>::same_point_test_sq(l0x, l0y, l1x, l1y));
  }

}


