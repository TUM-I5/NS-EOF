#include "StdAfx.hpp"

#include "MeshsizeFactory.hpp"

MeshsizeFactory& MeshsizeFactory::getInstance() {
  static MeshsizeFactory singleton;
  return singleton;
}

void MeshsizeFactory::initMeshsize(Parameters& parameters) {
  // Initialise meshsize
  switch (parameters.geometry.meshsizeType) {
  // Uniform meshsize
  case Uniform:
    parameters.meshsize = new UniformMeshsize(parameters);
    break;
  // Tanh-stretched mesh
  case TanhStretching:
    parameters.meshsize = new TanhMeshStretching(
      parameters,
      static_cast<bool>(parameters.geometry.stretchX),
      static_cast<bool>(parameters.geometry.stretchY),
      static_cast<bool>(parameters.geometry.stretchZ)
    );
    break;
  default:
    throw std::runtime_error("Unknown meshsize type!");
    break;
  }

  if (parameters.meshsize == NULL) {
    throw std::runtime_error("parameters.meshsize == NULL!");
  }
}
