#include "StdAfx.hpp"

#include "MovingWallStencils.hpp"

Stencils::MovingWallVelocityStencil::MovingWallVelocityStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters) {}

void Stencils::MovingWallVelocityStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0] = parameters_.walls.vectorLeft[0];
  flowField.getVelocity().getVector(i, j)[1] = 2 * parameters_.walls.vectorLeft[1]
                                               - flowField.getVelocity().getVector(i + 1, j)[1];
}

void Stencils::MovingWallVelocityStencil::applyRightWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i - 1, j)[0] = parameters_.walls.vectorRight[0];
  flowField.getVelocity().getVector(i, j)[1]     = 2 * parameters_.walls.vectorRight[1]
                                               - flowField.getVelocity().getVector(i - 1, j)[1];
}

void Stencils::MovingWallVelocityStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0] = 2 * parameters_.walls.vectorBottom[0]
                                               - flowField.getVelocity().getVector(i, j + 1)[0];
  flowField.getVelocity().getVector(i, j)[1] = parameters_.walls.vectorBottom[1];
}

void Stencils::MovingWallVelocityStencil::applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0] = 2 * parameters_.walls.vectorTop[0]
                                               - flowField.getVelocity().getVector(i, j - 1)[0];
  flowField.getVelocity().getVector(i, j - 1)[1] = parameters_.walls.vectorTop[1];
}

void Stencils::MovingWallVelocityStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = parameters_.walls.vectorLeft[0];
  flowField.getVelocity().getVector(i, j, k)[1] = 2 * parameters_.walls.vectorLeft[1]
                                                  - flowField.getVelocity().getVector(i + 1, j, k)[1];
  flowField.getVelocity().getVector(i, j, k)[2] = 2 * parameters_.walls.vectorLeft[2]
                                                  - flowField.getVelocity().getVector(i + 1, j, k)[2];
}

void Stencils::MovingWallVelocityStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i - 1, j, k)[0] = parameters_.walls.vectorRight[0];
  flowField.getVelocity().getVector(i, j, k)[1]     = 2 * parameters_.walls.vectorRight[1]
                                                  - flowField.getVelocity().getVector(i - 1, j, k)[1];
  flowField.getVelocity().getVector(i, j, k)[2] = 2 * parameters_.walls.vectorRight[2]
                                                  - flowField.getVelocity().getVector(i - 1, j, k)[2];
}

void Stencils::MovingWallVelocityStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = 2 * parameters_.walls.vectorBottom[0]
                                                  - flowField.getVelocity().getVector(i, j + 1, k)[0];
  flowField.getVelocity().getVector(i, j, k)[1] = parameters_.walls.vectorBottom[1];
  flowField.getVelocity().getVector(i, j, k)[2] = 2 * parameters_.walls.vectorBottom[2]
                                                  - flowField.getVelocity().getVector(i, j + 1, k)[2];
}

void Stencils::MovingWallVelocityStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = 2 * parameters_.walls.vectorTop[0]
                                                  - flowField.getVelocity().getVector(i, j - 1, k)[0];
  flowField.getVelocity().getVector(i, j - 1, k)[1] = parameters_.walls.vectorTop[1];
  flowField.getVelocity().getVector(i, j, k)[2]     = 2 * parameters_.walls.vectorTop[2]
                                                  - flowField.getVelocity().getVector(i, j - 1, k)[2];
}

void Stencils::MovingWallVelocityStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = 2 * parameters_.walls.vectorFront[0]
                                                  - flowField.getVelocity().getVector(i, j, k + 1)[0];
  flowField.getVelocity().getVector(i, j, k)[1] = 2 * parameters_.walls.vectorFront[1]
                                                  - flowField.getVelocity().getVector(i, j, k + 1)[1];
  flowField.getVelocity().getVector(i, j, k)[2] = parameters_.walls.vectorFront[2];
}

void Stencils::MovingWallVelocityStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = 2 * parameters_.walls.vectorBack[0]
                                                  - flowField.getVelocity().getVector(i, j, k - 1)[0];
  flowField.getVelocity().getVector(i, j, k)[1] = 2 * parameters_.walls.vectorBack[1]
                                                  - flowField.getVelocity().getVector(i, j, k - 1)[1];
  flowField.getVelocity().getVector(i, j, k - 1)[2] = parameters_.walls.vectorBack[2];
}

Stencils::MovingWallFGHStencil::MovingWallFGHStencil(const Parameters& parameters):
  BoundaryStencil<FlowField>(parameters) {}

void Stencils::MovingWallFGHStencil::applyLeftWall(FlowField& flowField, int i, int j) {
  flowField.getFGH().getVector(i, j)[0] = parameters_.walls.vectorLeft[0];
}

void Stencils::MovingWallFGHStencil::applyRightWall(FlowField& flowField, int i, int j) {
  flowField.getFGH().getVector(i - 1, j)[0] = parameters_.walls.vectorRight[0];
}

void Stencils::MovingWallFGHStencil::applyBottomWall(FlowField& flowField, int i, int j) {
  flowField.getFGH().getVector(i, j)[1] = parameters_.walls.vectorBottom[1];
}

void Stencils::MovingWallFGHStencil::applyTopWall(FlowField& flowField, int i, int j) {
  flowField.getFGH().getVector(i, j - 1)[1] = parameters_.walls.vectorTop[1];
}

void Stencils::MovingWallFGHStencil::applyLeftWall(FlowField& flowField, int i, int j, int k) {
  flowField.getFGH().getVector(i, j, k)[0] = parameters_.walls.vectorLeft[0];
}

void Stencils::MovingWallFGHStencil::applyRightWall(FlowField& flowField, int i, int j, int k) {
  flowField.getFGH().getVector(i - 1, j, k)[0] = parameters_.walls.vectorRight[0];
}

void Stencils::MovingWallFGHStencil::applyBottomWall(FlowField& flowField, int i, int j, int k) {
  flowField.getFGH().getVector(i, j, k)[1] = parameters_.walls.vectorBottom[1];
}

void Stencils::MovingWallFGHStencil::applyTopWall(FlowField& flowField, int i, int j, int k) {
  flowField.getFGH().getVector(i, j - 1, k)[1] = parameters_.walls.vectorTop[1];
}

void Stencils::MovingWallFGHStencil::applyFrontWall(FlowField& flowField, int i, int j, int k) {
  flowField.getFGH().getVector(i, j, k)[2] = parameters_.walls.vectorFront[2];
}

void Stencils::MovingWallFGHStencil::applyBackWall(FlowField& flowField, int i, int j, int k) {
  flowField.getFGH().getVector(i, j, k - 1)[2] = parameters_.walls.vectorBack[2];
}
