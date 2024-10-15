#include "StdAfx.hpp"

#include "LinearSolver.hpp"

Solvers::LinearSolver::LinearSolver(FlowField& flowField, const Parameters& parameters):
  flowField_(flowField),
  parameters_(parameters) {}
