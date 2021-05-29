#include <catch2/catch.hpp>

#include <limits>

#include "../src/geom2d.hpp"

using linegeom=geomalgos2d::linegeo2d<double, double>;
using namespace geomalgos2d;

TEST_CASE("closest point to non-degenerated line", "[closest][line]") {
  SECTION("sample out of line, closest between start/end") {
    double px = 0.0, py = 1.0;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::line);
    REQUIRE(hasClosest);
    REQUIRE(param == 0.5);
  }
  SECTION("sample out of line, closest before start") {
    double px = -1.0, py = 0;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::line);
    REQUIRE(hasClosest);
    REQUIRE(param == -0.5);
  }
  SECTION("sample out of line, closest after end") {
    double px = 0, py = 4;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::line);
    REQUIRE(hasClosest);
    REQUIRE(param == 2);
  }
  SECTION("sample on the line, between but not on start/end") {
    double px = 0.5, py = 0.5;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::line);
    REQUIRE(hasClosest);
    REQUIRE(param == 0.5);
  }
  SECTION("sample is line start") {
    double px = 0, py = 0;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::line);
    REQUIRE(hasClosest);
    REQUIRE(param == 0.0);
  }
  SECTION("sample is line end") {
    double px = 1, py = 1;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::line);
    REQUIRE(hasClosest);
    REQUIRE(param == 1.0);
  }
  SECTION("sample on the line, before start") {
    double px = -0.5, py = -0.5;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::line);
    REQUIRE(hasClosest);
    REQUIRE(param == -0.5);
  }
  SECTION("sample on the line, after end") {
    double px = 2, py = 2;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::line);
    REQUIRE(hasClosest);
    REQUIRE(param == 2);
  }
}

TEST_CASE("closest point to non-degenerated ray", "[closest][ray]") {
  SECTION("sample out of line, closest between start/end") {
    double px = 0.0, py = 1.0;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::ray);
    REQUIRE(hasClosest);
    REQUIRE(param == 0.5);
  }
  SECTION("sample out of line, closest before start") {
    double px = -1.0, py = 0;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::ray);
    REQUIRE(hasClosest);
    REQUIRE(param == 0); // closest is ray start
  }
  SECTION("sample out of line, closest after end") {
    double px = 0, py = 4;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::ray);
    REQUIRE(hasClosest);
    REQUIRE(param == 2);
  }
  SECTION("sample on the line, between but not on start/end") {
    double px = 0.5, py = 0.5;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::ray);
    REQUIRE(hasClosest);
    REQUIRE(param == 0.5);
  }
  SECTION("sample is line start") {
    double px = 0, py = 0;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::ray);
    REQUIRE(hasClosest);
    REQUIRE(param == 0.0);
  }
  SECTION("sample is line end") {
    double px = 1, py = 1;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::ray);
    REQUIRE(hasClosest);
    REQUIRE(param == 1.0);
  }
  SECTION("sample on the line, before start") {
    double px = -0.5, py = -0.5;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::ray);
    REQUIRE(hasClosest);
    REQUIRE(param == 0); // clamped to no negatives
  }
  SECTION("sample on the line, after end") {
    double px = 2, py = 2;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::ray);
    REQUIRE(hasClosest);
    REQUIRE(param == 2);
  }
}

TEST_CASE("closest point to non-degenerated segment", "[closest][segment]") {
  SECTION("sample out of line, closest between start/end") {
    double px = 0.0, py = 1.0;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::segment);
    REQUIRE(hasClosest);
    REQUIRE(param == 0.5);
  }
  SECTION("sample out of line, closest before start") {
    double px = -1.0, py = 0;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::segment);
    REQUIRE(hasClosest);
    REQUIRE(param == 0); // closest is segment start
  }
  SECTION("sample out of line, closest after end") {
    double px = 0, py = 4;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::segment);
    REQUIRE(hasClosest);
    REQUIRE(param == 1); // closest is segment end
  }
  SECTION("sample on the line, between but not on start/end") {
    double px = 0.5, py = 0.5;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::segment);
    REQUIRE(hasClosest);
    REQUIRE(param == 0.5);
  }
  SECTION("sample is line start") {
    double px = 0, py = 0;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::segment);
    REQUIRE(hasClosest);
    REQUIRE(param == 0.0);
  }
  SECTION("sample is line end") {
    double px = 1, py = 1;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::segment);
    REQUIRE(hasClosest);
    REQUIRE(param == 1.0);
  }
  SECTION("sample on the line, before start") {
    double px = -0.5, py = -0.5;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::segment);
    REQUIRE(hasClosest);
    REQUIRE(param == 0); // clamped to no negatives
  }
  SECTION("sample on the line, after end") {
    double px = 2, py = 2;
    double param = -1.0;
    bool hasClosest = linegeom::closest_point_param(px, py, 0.0, 0.0, 1.0, 1.0, param, linear_variety::segment);
    REQUIRE(hasClosest);
    REQUIRE(param == 1); //clamped to [0,1]
  }
}



