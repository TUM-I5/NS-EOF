#pragma once

#include "Definitions.hpp"
#include "Parameters.hpp"

/**
 * Initialises the meshsize in the Parameters. Must be called after configuring (Configuration and
 * PetscParallelConfiguration). We therefore make use of the singleton/factory pattern.
 */
class MeshsizeFactory {
private:
  MeshsizeFactory() = default;

public:
  ~MeshsizeFactory() = default;

  static MeshsizeFactory& getInstance();

  void initMeshsize(Parameters& parameters);
};
