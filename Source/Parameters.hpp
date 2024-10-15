#pragma once

#include "BoundaryType.hpp"
#include "Definitions.hpp"
#include "Meshsize.hpp"

//! Classes for the parts of the parameters
//@{
class TimestepParameters {
public:
  RealType dt  = 0; //! Timestep
  RealType tau = 0; //! Security factor
};

class SimulationParameters {
public:
  RealType    finalTime = 0; //! Final time for the simulation
  std::string type;          //! Type of the simulation (DNS vs. Turbulence)
  std::string scenario;      //! If channel or cavity, for example
};

class EnvironmentalParameters {
public:
  // Gravity components
  RealType gx = 0;
  RealType gy = 0;
  RealType gz = 0;
};

class FlowParameters {
public:
  RealType Re = 0; //! Reynolds number
};

class SolverParameters {
public:
  RealType gamma         = 0;  //! Donor cell balance coefficient
  int      maxIterations = -1; //! Maximum number of iterations in the linear solver
};

class GeometricParameters {
public:
  // Dimensions
  int dim = -1;

  // Number of cells
  int sizeX = -1;
  int sizeY = -1;
  int sizeZ = -1;

  // Cell sizing
  RealType lengthX = 0;
  RealType lengthY = 0;
  RealType lengthZ = 0;

  // Meshsize type
  int meshsizeType = -1;

  // For mesh-stretching
  int stretchX = -1;
  int stretchY = -1;
  int stretchZ = -1;
};

class WallParameters {
public:
  // Scalar value definition. Used to define the pressure, for example.
  RealType scalarLeft;
  RealType scalarRight;
  RealType scalarBottom;
  RealType scalarTop;
  RealType scalarFront;
  RealType scalarBack;

  // Vector values at the boundaries, to define, for example, the velocities.
  RealType vectorLeft[3];
  RealType vectorRight[3];
  RealType vectorBottom[3];
  RealType vectorTop[3];
  RealType vectorFront[3];
  RealType vectorBack[3];

  // Define how will the boundary behave
  BoundaryType typeLeft;
  BoundaryType typeRight;
  BoundaryType typeTop;
  BoundaryType typeBottom;
  BoundaryType typeFront;
  BoundaryType typeBack;
};

class VTKParameters {
public:
  RealType    interval = 0; //! Time interval for file printing
  std::string prefix;       //! Output filename
};

class StdOutParameters {
public:
  RealType interval = 0;
};

class ParallelParameters {
public:
  int rank = -1; //! Rank of the current processor

  int numProcessors[3]; //! Array with the number of processors in each direction

  //@brief Ranks of the neighbours
  //@{
  int leftNb   = MPI_PROC_NULL;
  int rightNb  = MPI_PROC_NULL;
  int bottomNb = MPI_PROC_NULL;
  int topNb    = MPI_PROC_NULL;
  int frontNb  = MPI_PROC_NULL;
  int backNb   = MPI_PROC_NULL;
  //@}

  int indices[3];     //! 3D indices to locate the array
  int localSize[3];   //! Size for the local flow field
  int firstCorner[3]; //! Position of the first element. Used for plotting.

#ifdef ENABLE_PETSC
  PetscInt* sizes[3]; //! Arrays with the sizes of the blocks in each direction.
#else
  int* sizes[3];
#endif
};

class BFStepParameters {
public:
  RealType xRatio = 0;
  RealType yRatio = 0;
};

//@}

/** A class to store and pass around the parameters
 */
class Parameters {
public:
  Parameters();
  ~Parameters();

  SimulationParameters    simulation;
  TimestepParameters      timestep;
  EnvironmentalParameters environment;
  FlowParameters          flow;
  SolverParameters        solver;
  GeometricParameters     geometry;
  WallParameters          walls;
  VTKParameters           vtk;
  ParallelParameters      parallel;
  StdOutParameters        stdOut;
  BFStepParameters        bfStep;
  // TODO WS2: include parameters for turbulence
  Meshsize* meshsize;
};
