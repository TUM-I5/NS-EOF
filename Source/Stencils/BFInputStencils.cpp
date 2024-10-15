#include "StdAfx.hpp"

#include "BFInputStencils.hpp"

RealType computeVelocity3D(int i, int j, int k, RealType stepSize, const Parameters& parameters) {
  const RealType posY = parameters.meshsize->getPosY(i, j, k);
  const RealType posZ = parameters.meshsize->getPosZ(i, j, k);
  const RealType dy   = parameters.meshsize->getDy(i, j, k);
  const RealType dz   = parameters.meshsize->getDz(i, j, k);

  if (posY + 0.5 * dy >= stepSize) {
    // Get the size of the inlet in Y. A 3 is subtracted because of the boundary cells.
    const RealType inletYSize = parameters.geometry.lengthY - stepSize;
    const RealType inletZSize = parameters.geometry.lengthZ;

    const RealType y = posY + 0.5 * dy - stepSize;
    const RealType z = posZ + 0.5 * dz;

    return 36.0 * parameters.walls.vectorLeft[0] / (inletZSize * inletZSize * inletYSize * inletYSize) * y
           * (y - inletYSize) * z * (z - inletZSize);
  } else {
    return 0.0;
  }
}

RealType computeVelocity2D(int i, int j, RealType stepSize, const Parameters& parameters) {
  const RealType posY = parameters.meshsize->getPosY(i, j);
  const RealType dy   = parameters.meshsize->getDy(i, j);

  if (posY + 0.5 * dy >= stepSize) {
    // Get the size of the inlet in Y. A 3 is subtracted because of the boundary cells.
    [[maybe_unused]] const RealType inletYSize = parameters.geometry.lengthY - stepSize;
    [[maybe_unused]] const RealType y = posY + 0.5 * dy - stepSize;

    return parameters.walls.vectorLeft[0]; // Uniform inflow
    //return 6.0 * parameters.walls.vectorLeft[0] / (inletYSize * inletYSize) * y * (inletYSize - y); // Parabolic inflow
  } else {
    return 0.0;
  }
}

Stencils::BFInputVelocityStencil::BFInputVelocityStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters)
  // Here, the obstacle size is set to zero if it was set as negative at the configuration.
  ,
  stepSize_(parameters.bfStep.yRatio > 0.0 ? parameters.bfStep.yRatio * parameters.geometry.lengthY : 0.0) {

  if (parameters_.geometry.dim == 2) {
    RealType posY   = parameters_.meshsize->getPosY(0, 0);
    RealType dy     = parameters_.meshsize->getDy(0, 0);
    RealType nextDy = parameters_.meshsize->getDy(0, 1);

    for (int j = 0; j < parameters_.geometry.sizeY - 1; ++j) {
      posY   = parameters_.meshsize->getPosY(0, j);
      dy     = parameters_.meshsize->getDy(0, j);
      nextDy = parameters_.meshsize->getDy(0, j + 1);

      // Check if stepSize is in this cell
      if (posY + 0.5 * dy < stepSize_ && stepSize_ <= posY + dy + 0.5 * nextDy) {
        stepSize_ = posY + dy;
        break;
      }
    }
  } else if (parameters_.geometry.dim == 3) {
    RealType posY   = parameters_.meshsize->getPosY(0, 0, 0);
    RealType dy     = parameters_.meshsize->getDy(0, 0, 0);
    RealType nextDy = parameters_.meshsize->getDy(0, 1, 0);

    for (int j = 0; j < parameters_.geometry.sizeY - 1; ++j) {
      posY   = parameters_.meshsize->getPosY(0, j, 0);
      dy     = parameters_.meshsize->getDy(0, j, 0);
      nextDy = parameters_.meshsize->getDy(0, j + 1, 0);

      if (posY + 0.5 * dy < stepSize_ && stepSize_ <= posY + dy + 0.5 * nextDy) {
        stepSize_ = posY + dy;
        break;
      }
    }
  }
}

void Stencils::BFInputVelocityStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0] = computeVelocity2D(i, j, stepSize_, parameters_);
  flowField.getVelocity().getVector(i, j)[1] = -flowField.getVelocity().getVector(i + 1, j)[1];
}

// Most of the functions are empty, and they shouldn't be called, assuming that the input is always located at the left.
void Stencils::BFInputVelocityStencil::applyRightWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}
void Stencils::BFInputVelocityStencil::applyBottomWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}
void Stencils::BFInputVelocityStencil::applyTopWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}

void Stencils::BFInputVelocityStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = computeVelocity3D(i, j, k, stepSize_, parameters_);
  flowField.getVelocity().getVector(i, j, k)[1] = -flowField.getVelocity().getVector(i + 1, j, k)[1];
  flowField.getVelocity().getVector(i, j, k)[2] = -flowField.getVelocity().getVector(i + 1, j, k)[2];
}

void Stencils::BFInputVelocityStencil::applyRightWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputVelocityStencil::applyBottomWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputVelocityStencil::applyTopWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputVelocityStencil::applyFrontWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
void Stencils::BFInputVelocityStencil::applyBackWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}

Stencils::BFInputFGHStencil::BFInputFGHStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters),
  stepSize_(parameters.bfStep.yRatio > 0.0 ? parameters.bfStep.yRatio * parameters.geometry.lengthY : 0.0) {}

void Stencils::BFInputFGHStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getFGH().getVector(i, j)[0] = computeVelocity2D(i, j, stepSize_, parameters_);
}

void Stencils::BFInputFGHStencil::applyRightWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}
void Stencils::BFInputFGHStencil::applyBottomWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}
void Stencils::BFInputFGHStencil::applyTopWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j
) {}

void Stencils::BFInputFGHStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getFGH().getVector(i, j, k)[0] = computeVelocity3D(i, j, k, stepSize_, parameters_);
}

void Stencils::BFInputFGHStencil::applyRightWall(
  [[maybe_unused]] FlowField&           flowField,
  [[maybe_unused]] int                  i,
  [[maybe_unused]] int                  j,
  [[maybe_unused]] int                  k
) {}
void Stencils::BFInputFGHStencil::applyBottomWall(
  [[maybe_unused]] FlowField&           flowField,
  [[maybe_unused]] int                  i,
  [[maybe_unused]] int                  j,
  [[maybe_unused]] int                  k
) {}
void Stencils::BFInputFGHStencil::applyTopWall(
  [[maybe_unused]] FlowField&           flowField,
  [[maybe_unused]] int                  i,
  [[maybe_unused]] int                  j,
  [[maybe_unused]] int                  k
) {}
void Stencils::BFInputFGHStencil::applyFrontWall(
  [[maybe_unused]] FlowField&           flowField,
  [[maybe_unused]] int                  i,
  [[maybe_unused]] int                  j,
  [[maybe_unused]] int                  k
) {}
void Stencils::BFInputFGHStencil::applyBackWall(
  [[maybe_unused]] FlowField& flowField, [[maybe_unused]] int i, [[maybe_unused]] int j, [[maybe_unused]] int k
) {}
