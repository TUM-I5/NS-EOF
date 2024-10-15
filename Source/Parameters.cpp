#include "StdAfx.hpp"

#include "Parameters.hpp"

Parameters::Parameters():
  simulation{},
  timestep{},
  environment{},
  flow{},
  solver{},
  geometry{},
  walls{},
  vtk{},
  parallel{},
  stdOut{},
  bfStep{},
  meshsize(NULL) {
}

Parameters::~Parameters() {
  if (meshsize != NULL) {
    delete meshsize;
    meshsize = NULL;
  }
}
