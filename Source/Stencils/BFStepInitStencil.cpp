#include "StdAfx.hpp"

#include "BFStepInitStencil.hpp"

Stencils::BFStepInitStencil::BFStepInitStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters),
  xLimit_(parameters.bfStep.xRatio * parameters.geometry.lengthX),
  yLimit_(parameters.bfStep.yRatio * parameters.geometry.lengthY) {}

void Stencils::BFStepInitStencil::apply(FlowField& flowField, int i, int j) {
  IntScalarField& flags  = flowField.getFlags();
  const RealType  posX   = parameters_.meshsize->getPosX(i, j);
  const RealType  posY   = parameters_.meshsize->getPosY(i, j);
  const RealType  dx     = parameters_.meshsize->getDx(i, j);
  const RealType  dy     = parameters_.meshsize->getDy(i, j);
  const RealType  nextDx = parameters_.meshsize->getDx(i + 1, j);
  const RealType  nextDy = parameters_.meshsize->getDy(i, j + 1);
  const RealType  lastDx = parameters_.meshsize->getDx(i - 1, j);
  const RealType  lastDy = parameters_.meshsize->getDy(i, j - 1);

  if (posX + 0.5 * dx < xLimit_ && posY + 0.5 * dy < yLimit_) {
    flags.getValue(i, j) = OBSTACLE_SELF;
  }
  if (posX - 0.5 * lastDx < xLimit_ && posY + 0.5 * dy < yLimit_) {
    flags.getValue(i, j) += OBSTACLE_LEFT;
  }
  if (posX + dx + 0.5 * nextDx < xLimit_ && posY + 0.5 * dy < yLimit_) {
    flags.getValue(i, j) += OBSTACLE_RIGHT;
  }
  if (posX + 0.5 * dx < xLimit_ && posY - 0.5 * lastDy < yLimit_) {
    flags.getValue(i, j) += OBSTACLE_BOTTOM;
  }
  if (posX + 0.5 * dx < xLimit_ && posY + dy + 0.5 * nextDy < yLimit_) {
    flags.getValue(i, j) += OBSTACLE_TOP;
  }
}

void Stencils::BFStepInitStencil::apply(FlowField& flowField, int i, int j, int k) {
  IntScalarField& flags  = flowField.getFlags();
  const RealType  posX   = parameters_.meshsize->getPosX(i, j, k);
  const RealType  posY   = parameters_.meshsize->getPosY(i, j, k);
  const RealType  dx     = parameters_.meshsize->getDx(i, j, k);
  const RealType  dy     = parameters_.meshsize->getDy(i, j, k);
  const RealType  nextDx = parameters_.meshsize->getDx(i + 1, j, k);
  const RealType  nextDy = parameters_.meshsize->getDy(i, j + 1, k);
  const RealType  lastDx = parameters_.meshsize->getDx(i - 1, j, k);
  const RealType  lastDy = parameters_.meshsize->getDy(i, j - 1, k);

  if (posX + 0.5 * dx < xLimit_ && posY + 0.5 * dy < yLimit_) {
    flags.getValue(i, j, k) = OBSTACLE_SELF;
    // The obstacle is 2D, so we can say the following:
    flags.getValue(i, j, k) += OBSTACLE_FRONT;
    flags.getValue(i, j, k) += OBSTACLE_BACK;
  }
  if (posX - 0.5 * lastDx < xLimit_ && posY + 0.5 * dy < yLimit_) {
    flags.getValue(i, j, k) += OBSTACLE_LEFT;
  }
  if (posX + dx + 0.5 * nextDx < xLimit_ && posY + 0.5 * dy < yLimit_) {
    flags.getValue(i, j, k) += OBSTACLE_RIGHT;
  }
  if (posX + 0.5 * dx < xLimit_ && posY - 0.5 * lastDy < yLimit_) {
    flags.getValue(i, j, k) += OBSTACLE_BOTTOM;
  }
  if (posX + 0.5 * dx < xLimit_ && posY + dx + 0.5 * nextDy < yLimit_) {
    flags.getValue(i, j, k) += OBSTACLE_TOP;
  }
}
