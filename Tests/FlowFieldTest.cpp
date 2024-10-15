#include "StdAfx.hpp"

#include <catch2/catch_test_macros.hpp>

#include "DataStructures.hpp"
#include "FlowField.hpp"

constexpr auto SIZE_X = 20;
constexpr auto SIZE_Y = 25;

TEST_CASE("Test flow field", "[single-file]") {
  spdlog::info("Testing flow field");

  FlowField field(SIZE_X, SIZE_Y);

  REQUIRE(field.getPressure().getScalar(10, 10) == 0);

  field.getPressure().getScalar(10, 10) = 10;
  REQUIRE(field.getPressure().getScalar(10, 10) / 2 == 5);

  REQUIRE(field.getFlags().getNx() == 23);

  REQUIRE(field.getFlags().getNy() == 28);

  field.getFlags().getValue(10, 10) = 7;

  REQUIRE(field.getFlags().getValue(10, 10) / 3 == 2);

  spdlog::info("Test for flow field completed successfully");
}
