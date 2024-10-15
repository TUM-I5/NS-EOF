#include "StdAfx.hpp"

#include "MaxUStencil.hpp"

Stencils::MaxUStencil::MaxUStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters),
  BoundaryStencil<FlowField>(parameters) {

  reset();
}

void Stencils::MaxUStencil::apply(FlowField& flowField, int i, int j) { cellMaxValue(flowField, i, j); }

void Stencils::MaxUStencil::apply(FlowField& flowField, int i, int j, int k) { cellMaxValue(flowField, i, j, k); }

void Stencils::MaxUStencil::applyLeftWall(FlowField& flowField, int i, int j) { cellMaxValue(flowField, i, j); }

void Stencils::MaxUStencil::applyRightWall(FlowField& flowField, int i, int j) { cellMaxValue(flowField, i, j); }

void Stencils::MaxUStencil::applyBottomWall(FlowField& flowField, int i, int j) { cellMaxValue(flowField, i, j); }

void Stencils::MaxUStencil::applyTopWall(FlowField& flowField, int i, int j) { cellMaxValue(flowField, i, j); }

void Stencils::MaxUStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  cellMaxValue(flowField, i, j, k);
}

void Stencils::MaxUStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  cellMaxValue(flowField, i, j, k);
}

void Stencils::MaxUStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  cellMaxValue(flowField, i, j, k);
}

void Stencils::MaxUStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  cellMaxValue(flowField, i, j, k);
}

void Stencils::MaxUStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  cellMaxValue(flowField, i, j, k);
}

void Stencils::MaxUStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  cellMaxValue(flowField, i, j, k);
}

void Stencils::MaxUStencil::cellMaxValue(FlowField& flowField, int i, int j) {
  RealType*      velocity = flowField.getVelocity().getVector(i, j);
  const RealType dx       = FieldStencil<FlowField>::parameters_.meshsize->getDx(i, j);
  const RealType dy       = FieldStencil<FlowField>::parameters_.meshsize->getDy(i, j);
  if (fabs(velocity[0]) / dx > maxValues_[0]) {
    maxValues_[0] = fabs(velocity[0]) / dx;
  }
  if (fabs(velocity[1]) / dy > maxValues_[1]) {
    maxValues_[1] = fabs(velocity[1]) / dy;
  }
}

void Stencils::MaxUStencil::cellMaxValue(FlowField& flowField, int i, int j, int k) {
  RealType*      velocity = flowField.getVelocity().getVector(i, j, k);
  const RealType dx       = FieldStencil<FlowField>::parameters_.meshsize->getDx(i, j, k);
  const RealType dy       = FieldStencil<FlowField>::parameters_.meshsize->getDy(i, j, k);
  const RealType dz       = FieldStencil<FlowField>::parameters_.meshsize->getDz(i, j, k);
  if (fabs(velocity[0]) / dx > maxValues_[0]) {
    maxValues_[0] = fabs(velocity[0]) / dx;
  }
  if (fabs(velocity[1]) / dy > maxValues_[1]) {
    maxValues_[1] = fabs(velocity[1]) / dy;
  }
  if (fabs(velocity[2]) / dz > maxValues_[2]) {
    maxValues_[2] = fabs(velocity[2]) / dz;
  }
}

void Stencils::MaxUStencil::reset() {
  maxValues_[0] = 0;
  maxValues_[1] = 0;
  maxValues_[2] = 0;
}

const RealType* Stencils::MaxUStencil::getMaxValues() const { return maxValues_; }
