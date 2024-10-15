#pragma once

#include "Definitions.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Solvers {

  // Abstract class for linear solvers for the pressure
  class LinearSolver {
  protected:
    FlowField&        flowField_;
    const Parameters& parameters_;

  public:
    LinearSolver(FlowField& flowField, const Parameters& parameters);
    virtual ~LinearSolver() = default;

    virtual void        solve() = 0;
    virtual inline void reInitMatrix() {}
  };

} // namespace Solvers
