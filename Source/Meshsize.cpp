#include "StdAfx.hpp"

#include "Meshsize.hpp"

#include "Parameters.hpp"

UniformMeshsize::UniformMeshsize(const Parameters& parameters):
  Meshsize(),
  dx_(parameters.geometry.lengthX / parameters.geometry.sizeX),
  dy_(parameters.geometry.lengthY / parameters.geometry.sizeY),
  dz_(parameters.geometry.dim == 3 ? parameters.geometry.lengthZ / parameters.geometry.sizeZ : 0.0),
  firstCornerX_(parameters.parallel.firstCorner[0]),
  firstCornerY_(parameters.parallel.firstCorner[1]),
  firstCornerZ_(parameters.geometry.dim == 3 ? parameters.parallel.firstCorner[2] : 0) {

  if (dx_ <= 0.0) {
    throw std::runtime_error("dx <= 0.0!");
  }
  if (dy_ <= 0.0) {
    throw std::runtime_error("dy <= 0.0!");
  }
  if (parameters.geometry.dim == 3) {
    if (dz_ <= 0.0) {
      throw std::runtime_error("dz <= 0.0!");
    }
  }
}

TanhMeshStretching::TanhMeshStretching(const Parameters& parameters, bool stretchX, bool stretchY, bool stretchZ):
  Meshsize(),
  uniformMeshsize_(parameters),
  lengthX_(parameters.geometry.lengthX),
  lengthY_(parameters.geometry.lengthY),
  lengthZ_(parameters.geometry.dim == 3 ? parameters.geometry.lengthZ : 0.0),
  sizeX_(parameters.geometry.sizeX),
  sizeY_(parameters.geometry.sizeY),
  sizeZ_(parameters.geometry.dim == 3 ? parameters.geometry.sizeZ : 1),
  firstCornerX_(parameters.parallel.firstCorner[0]),
  firstCornerY_(parameters.parallel.firstCorner[1]),
  firstCornerZ_(parameters.geometry.dim == 3 ? parameters.parallel.firstCorner[2] : 0),
  stretchX_(stretchX),
  stretchY_(stretchY),
  stretchZ_(stretchZ),
  deltaS_(2.7),
  tanhDeltaS_(tanh(2.7)) // This parameters is chosen as 2.7 as used also in the dissertation by Tobias Neckel
  ,
  dxMin_(
    stretchX
      ? 0.5 * parameters.geometry.lengthX * (1.0 + tanh(deltaS_ * (2.0 / sizeX_ - 1.0)) / tanhDeltaS_)
      : uniformMeshsize_.getDx(0, 0)
  ),
  dyMin_(
    stretchY
      ? 0.5 * parameters.geometry.lengthY * (1.0 + tanh(deltaS_ * (2.0 / sizeY_ - 1.0)) / tanhDeltaS_)
      : uniformMeshsize_.getDy(0, 0)
  ),
  dzMin_(
    stretchZ
      ? 0.5 * parameters.geometry.lengthZ * (1.0 + tanh(deltaS_ * (2.0 / sizeZ_ - 1.0)) / tanhDeltaS_)
      : uniformMeshsize_.getDz(0, 0, 0)
  ) {}
