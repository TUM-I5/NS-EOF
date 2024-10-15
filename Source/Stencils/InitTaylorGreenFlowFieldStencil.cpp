#include "StdAfx.hpp"

#include "InitTaylorGreenFlowFieldStencil.hpp"

Stencils::InitTaylorGreenFlowFieldStencil::InitTaylorGreenFlowFieldStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters),
  pi2_(2.0 * 3.141592653589793238),
  domainSize_(initializeDomainSize(parameters)) {}

Stencils::InitTaylorGreenFlowFieldStencil::~InitTaylorGreenFlowFieldStencil() {
  if (domainSize_ != NULL) {
    delete[] domainSize_;
  }
}

void Stencils::InitTaylorGreenFlowFieldStencil::apply(FlowField& flowField, int i, int j) {
  RealType        coords[3] = {0.0, 0.0, 0.0};
  RealType* const velocity  = flowField.getVelocity().getVector(i, j);
  computeGlobalCoordinates(coords, i, j);
  // Initialize velocities
  velocity[0] = sin(pi2_ * (coords[0] + 0.5 * parameters_.meshsize->getDx(i, j)) / domainSize_[0])
                * sin(pi2_ * coords[1] / domainSize_[1]);
  velocity[1] = cos(pi2_ * coords[0] / domainSize_[0])
                * cos(pi2_ * (coords[1] + 0.5 * parameters_.meshsize->getDy(i, j)) / domainSize_[1]);
}

void Stencils::InitTaylorGreenFlowFieldStencil::apply(FlowField& flowField, int i, int j, int k) {
  RealType        coords[3] = {0.0, 0.0, 0.0};
  RealType* const velocity  = flowField.getVelocity().getVector(i, j, k);
  computeGlobalCoordinates(coords, i, j, k);
  // Initialize velocities
  velocity[0] = cos(pi2_ * (coords[0] + 0.5 * parameters_.meshsize->getDx(i, j, k)) / domainSize_[0])
                * sin(pi2_ * coords[1] / domainSize_[1]) * sin(pi2_ * coords[2] / domainSize_[2]);
  velocity[1] = sin(pi2_ * coords[0] / domainSize_[0])
                * cos(pi2_ * (coords[1] + 0.5 * parameters_.meshsize->getDy(i, j, k)) / domainSize_[1])
                * sin(pi2_ * coords[2] / domainSize_[2]);
  velocity[2] = sin(pi2_ * coords[0] / domainSize_[0]) * sin(pi2_ * coords[1] / domainSize_[1])
                * cos(pi2_ * (coords[2] + 0.5 * parameters_.meshsize->getDz(i, j, k)) / domainSize_[2]);
}

RealType* Stencils::InitTaylorGreenFlowFieldStencil::initializeDomainSize(const Parameters& parameters) const {
  RealType* domainSize = new RealType[3];
  if (domainSize == NULL) {
    throw std::runtime_error("domainSize == NULL");
  }
  domainSize[0] = parameters.geometry.lengthX;
  domainSize[1] = parameters.geometry.lengthY;
  domainSize[2] = parameters.geometry.lengthZ;
  return domainSize;
}

void Stencils::InitTaylorGreenFlowFieldStencil::computeGlobalCoordinates(RealType* coords, int i, int j, int k) const {
  coords[0] = parameters_.meshsize->getPosX(i, j, k) + 0.5 * parameters_.meshsize->getDx(i, j, k);
  coords[1] = parameters_.meshsize->getPosY(i, j, k) + 0.5 * parameters_.meshsize->getDy(i, j, k);
  if (parameters_.geometry.dim == 3) {
    coords[2] = parameters_.meshsize->getPosZ(i, j, k) + 0.5 * parameters_.meshsize->getDz(i, j, k);
  }
}
