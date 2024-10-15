#pragma once

#include "LinearSolver.hpp"

namespace Solvers {

  class SORSolver: public LinearSolver {
  public:
    SORSolver(FlowField& flowField, const Parameters& parameters);
    ~SORSolver() override = default;

    void solve() override;
  };

} // namespace Solvers
