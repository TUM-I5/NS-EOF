#pragma once

// Types of boundary
enum BoundaryType {
  DIRICHLET, // Also used for the case of non-moving wall
  PERIODIC,
  PARALLEL_BOUNDARY,
  NEUMANN
};

static constexpr int OBSTACLE_SELF   = 1 << 0;
static constexpr int OBSTACLE_LEFT   = 1 << 1;
static constexpr int OBSTACLE_RIGHT  = 1 << 2;
static constexpr int OBSTACLE_BOTTOM = 1 << 3;
static constexpr int OBSTACLE_TOP    = 1 << 4;
static constexpr int OBSTACLE_FRONT  = 1 << 5;
static constexpr int OBSTACLE_BACK   = 1 << 6;
