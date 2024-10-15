#include "StdAfx.hpp"

#include "RHSStencil.hpp"

Stencils::RHSStencil::RHSStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters) {}

void Stencils::RHSStencil::apply(FlowField& flowField, int i, int j) {
  flowField.getRHS().getScalar(i, j) = 1.0 / parameters_.timestep.dt *
        ((flowField.getFGH().getVector(i, j)[0] - flowField.getFGH().getVector(i - 1, j)[0]) / parameters_.meshsize->getDx(i, j) +
         (flowField.getFGH().getVector(i, j)[1] - flowField.getFGH().getVector(i, j - 1)[1]) / parameters_.meshsize->getDy(i, j));
}

void Stencils::RHSStencil::apply(FlowField& flowField, int i, int j, int k) {
  flowField.getRHS().getScalar(i, j, k) = 1.0 / parameters_.timestep.dt *
        ((flowField.getFGH().getVector(i, j, k)[0] - flowField.getFGH().getVector(i - 1, j, k)[0]) / parameters_.meshsize->getDx(i, j, k) +
         (flowField.getFGH().getVector(i, j, k)[1] - flowField.getFGH().getVector(i, j - 1, k)[1]) / parameters_.meshsize->getDy(i, j, k) +
         (flowField.getFGH().getVector(i, j, k)[2] - flowField.getFGH().getVector(i, j, k - 1)[2]) / parameters_.meshsize->getDz(i, j, k));
}
