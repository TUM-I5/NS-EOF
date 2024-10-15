#include "StdAfx.hpp"

#include <catch2/catch_test_macros.hpp>

#include "DataStructures.hpp"

constexpr auto SIZE_X = 10;
constexpr auto SIZE_Y = 10;
constexpr auto SIZE_Z = 10;

TEST_CASE("Test scalar field", "[single-file]") {
  spdlog::info("Testing scalar fields");

  ScalarField sfield2D(SIZE_X, SIZE_Y);
  ScalarField sfield3D(SIZE_X, SIZE_Y, SIZE_Z);

  RealType d2counter = 1.0;
  RealType d3counter = 1.0;

  // Fill the fields completely with stuff
  for (int i = 0; i < SIZE_X; i++) {
    for (int j = 0; j < SIZE_Y; j++) {
      sfield2D.getScalar(i, j) = d2counter;
      d2counter += 0.5;

      for (int k = 0; k < SIZE_Z; k++) {
        sfield3D.getScalar(i, j, k) = d3counter;
        d3counter += 0.5;
      }
    }
  }

  // Now read and see if the same values come out
  d2counter = 1.0;
  d3counter = 1.0;

  for (int i = 0; i < SIZE_X; i++) {
    for (int j = 0; j < SIZE_Y; j++) {
      REQUIRE(sfield2D.getScalar(i, j) == d2counter);
      d2counter += 0.5;

      for (int k = 0; k < SIZE_Z; k++) {
        REQUIRE(sfield3D.getScalar(i, j, k) == d3counter);
        d3counter += 0.5;
      }
    }
  }

  spdlog::info("Test for scalar fields completed successfully");
}