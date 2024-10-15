#include "StdAfx.hpp"

#include "FlowField.hpp"

FlowField::FlowField(int Nx, int Ny):
  sizeX_(Nx),
  sizeY_(Ny),
  sizeZ_(1),
  cellsX_(Nx + 3),
  cellsY_(Ny + 3),
  cellsZ_(1)
  // Pressure field doesn't need to have an extra layer, but this allows to address the same
  // positions with the same iterator for both pressures and velocities.
  ,
  pressure_(ScalarField(Nx + 3, Ny + 3)),
  velocity_(VectorField(Nx + 3, Ny + 3)),
  flags_(IntScalarField(Nx + 3, Ny + 3)),
  FGH_(VectorField(Nx + 3, Ny + 3)),
  RHS_(ScalarField(Nx + 3, Ny + 3)) {

  ASSERTION(Nx > 0);
  ASSERTION(Ny > 0);
}

FlowField::FlowField(int Nx, int Ny, int Nz):
  sizeX_(Nx),
  sizeY_(Ny),
  sizeZ_(Nz),
  cellsX_(Nx + 3),
  cellsY_(Ny + 3),
  cellsZ_(Nz + 3),
  pressure_(ScalarField(Nx + 3, Ny + 3, Nz + 3)),
  velocity_(VectorField(Nx + 3, Ny + 3, Nz + 3)),
  flags_(IntScalarField(Nx + 3, Ny + 3, Nz + 3)),
  FGH_(VectorField(Nx + 3, Ny + 3, Nz + 3)),
  RHS_(ScalarField(Nx + 3, Ny + 3, Nz + 3)) {

  ASSERTION(Nx > 0);
  ASSERTION(Ny > 0);
  ASSERTION(Nz > 0);
}

FlowField::FlowField(const Parameters& parameters):
  sizeX_(parameters.parallel.localSize[0]),
  sizeY_(parameters.parallel.localSize[1]),
  sizeZ_(parameters.parallel.localSize[2]),
  cellsX_(sizeX_ + 3),
  cellsY_(sizeY_ + 3),
  cellsZ_(parameters.geometry.dim == 2 ? 1 : sizeZ_ + 3),
  pressure_(
    parameters.geometry.dim == 2 ? ScalarField(sizeX_ + 3, sizeY_ + 3) : ScalarField(sizeX_ + 3, sizeY_ + 3, sizeZ_ + 3)
  ),
  velocity_(
    parameters.geometry.dim == 2 ? VectorField(sizeX_ + 3, sizeY_ + 3) : VectorField(sizeX_ + 3, sizeY_ + 3, sizeZ_ + 3)
  ),
  flags_(
    parameters.geometry.dim == 2
      ? IntScalarField(sizeX_ + 3, sizeY_ + 3)
      : IntScalarField(sizeX_ + 3, sizeY_ + 3, sizeZ_ + 3)
  ),
  FGH_(
    parameters.geometry.dim == 2 ? VectorField(sizeX_ + 3, sizeY_ + 3) : VectorField(sizeX_ + 3, sizeY_ + 3, sizeZ_ + 3)
  ),
  RHS_(
    parameters.geometry.dim == 2 ? ScalarField(sizeX_ + 3, sizeY_ + 3) : ScalarField(sizeX_ + 3, sizeY_ + 3, sizeZ_ + 3)
  ) {}

int FlowField::getNx() const { return sizeX_; }

int FlowField::getNy() const { return sizeY_; }

int FlowField::getNz() const { return sizeZ_; }

int FlowField::getCellsX() const { return cellsX_; }

int FlowField::getCellsY() const { return cellsY_; }

int FlowField::getCellsZ() const { return cellsZ_; }

ScalarField& FlowField::getPressure() { return pressure_; }

VectorField& FlowField::getVelocity() { return velocity_; }

IntScalarField& FlowField::getFlags() { return flags_; }

VectorField& FlowField::getFGH() { return FGH_; }

ScalarField& FlowField::getRHS() { return RHS_; }

void FlowField::getPressureAndVelocity(RealType& pressure, RealType* const velocity, int i, int j) {
  RealType* vHere = getVelocity().getVector(i, j);
  RealType* vLeft = getVelocity().getVector(i - 1, j);
  RealType* vDown = getVelocity().getVector(i, j - 1);

  velocity[0] = (vHere[0] + vLeft[0]) / 2;
  velocity[1] = (vHere[1] + vDown[1]) / 2;

  pressure = getPressure().getScalar(i, j);
}

void FlowField::getPressureAndVelocity(RealType& pressure, RealType* const velocity, int i, int j, int k) {
  RealType* vHere = getVelocity().getVector(i, j, k);
  RealType* vLeft = getVelocity().getVector(i - 1, j, k);
  RealType* vDown = getVelocity().getVector(i, j - 1, k);
  RealType* vBack = getVelocity().getVector(i, j, k - 1);

  velocity[0] = (vHere[0] + vLeft[0]) / 2;
  velocity[1] = (vHere[1] + vDown[1]) / 2;
  velocity[2] = (vHere[2] + vBack[2]) / 2;

  pressure = getPressure().getScalar(i, j, k);
}
