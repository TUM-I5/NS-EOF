#include "StdAfx.hpp"

#ifdef ENABLE_PETSC

#include "PetscSolver.hpp"

static constexpr unsigned char LEFT_WALL_BIT   = 1 << 0;
static constexpr unsigned char RIGHT_WALL_BIT  = 1 << 1;
static constexpr unsigned char BOTTOM_WALL_BIT = 1 << 2;
static constexpr unsigned char TOP_WALL_BIT    = 1 << 3;
static constexpr unsigned char FRONT_WALL_BIT  = 1 << 4;
static constexpr unsigned char BACK_WALL_BIT   = 1 << 5;

// This function returns the ranges to work on the pressure with the non-boundary stencil.
// Since the domain PETSc deals with has an additional layer of cells, the size is clipped to
// ignore them.
void createLimits(Parameters& parameters, DM& da, int* limitsX, int* limitsY, int* limitsZ) {
  // Location of the first element and sizes of the subdomain in each dimension.
  PetscInt firstX, lengthX, firstY, lengthY, firstZ, lengthZ;

  DMDAGetCorners(da, &firstX, &firstY, &firstZ, &lengthX, &lengthY, &lengthZ);

  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // Preliminary set for the iteration domain
  limitsX[0] = firstX;
  limitsX[1] = firstX + lengthX;
  limitsY[0] = firstY;
  limitsY[1] = firstY + lengthY;
  limitsZ[0] = firstZ;
  limitsZ[1] = firstZ + lengthZ;

  // The indices are used to determine the local range of iterations. If the domain is not in a
  // global boundary, we don't have to consider special conditions for the boundaries. Otherwise,
  // these are shifted so that the body iterations don't touch the boundaries.

  // Check left wall
  if (parameters.parallel.indices[0] == 0) {
    limitsX[0]++;
  }

  // Check right wall
  if (parameters.parallel.indices[0] == parameters.parallel.numProcessors[0] - 1) {
    limitsX[1]--;
  }

  // Check bottom wall
  if (parameters.parallel.indices[1] == 0) {
    limitsY[0]++;
  }

  // Check top wall
  if (parameters.parallel.indices[1] == parameters.parallel.numProcessors[1] - 1) {
    limitsY[1]--;
  }

  // Check front wall
  if (parameters.parallel.indices[2] == 0) {
    limitsZ[0]++;
  }

  // Check back wall
  if (parameters.parallel.indices[2] == parameters.parallel.numProcessors[2] - 1) {
    limitsZ[1]--;
  }
}

Solvers::PetscUserCtx::PetscUserCtx(Parameters& parameters, FlowField& flowField):
  parameters_(parameters),
  flowField_(flowField) {}

Parameters& Solvers::PetscUserCtx::getParameters() { return parameters_; }

FlowField& Solvers::PetscUserCtx::getFlowField() { return flowField_; }

PetscErrorCode computeMatrix2D(KSP ksp, Mat A, Mat pc, void* ctx);
PetscErrorCode computeMatrix3D(KSP ksp, Mat A, Mat pc, void* ctx);

PetscErrorCode computeRHS2D(KSP ksp, Vec b, void* ctx);
PetscErrorCode computeRHS3D(KSP ksp, Vec b, void* ctx);

Solvers::PetscSolver::PetscSolver(FlowField& flowField, Parameters& parameters):
  LinearSolver(flowField, parameters),
  ctx_(parameters, flowField) {

  // Set the type of boundary nodes of the system
  DMBoundaryType bx = DM_BOUNDARY_NONE, by = DM_BOUNDARY_NONE, bz = DM_BOUNDARY_NONE;

  if (parameters.walls.typeLeft == PERIODIC) {
    bx = DM_BOUNDARY_PERIODIC;
  }
  if (parameters.walls.typeBottom == PERIODIC) {
    by = DM_BOUNDARY_PERIODIC;
  }
  if (parameters.walls.typeFront == PERIODIC) {
    bz = DM_BOUNDARY_PERIODIC;
  }

  KSPCreate(PETSC_COMM_WORLD, &ksp_);
  PCCreate(PETSC_COMM_WORLD, &pc_);

  PetscErrorCode (*computeMatrix)(KSP, Mat, Mat, void*) = NULL;

  if (parameters_.geometry.dim == 2) {
    computeMatrix = computeMatrix2D;
    DMDACreate2d(
      PETSC_COMM_WORLD,
      bx,
      by,
      DMDA_STENCIL_STAR,
      parameters_.geometry.sizeX + 2,
      parameters.geometry.sizeY + 2,
      parameters_.parallel.numProcessors[0],
      parameters_.parallel.numProcessors[1],
      1,
      2,
      parameters_.parallel.sizes[0],
      parameters_.parallel.sizes[1],
      &da_
    );
  } else if (parameters_.geometry.dim == 3) {
    computeMatrix = computeMatrix3D;
    DMDACreate3d(
      PETSC_COMM_WORLD,
      bx,
      by,
      bz,
      DMDA_STENCIL_STAR,
      parameters_.geometry.sizeX + 2,
      parameters.geometry.sizeY + 2,
      parameters_.geometry.sizeZ + 2,
      parameters_.parallel.numProcessors[0],
      parameters_.parallel.numProcessors[1],
      parameters_.parallel.numProcessors[2],
      1,
      2, // Degrees of freedom and stencil length
      parameters_.parallel.sizes[0],
      parameters_.parallel.sizes[1],
      parameters_.parallel.sizes[2],
      &da_
    );
  }

  DMSetFromOptions(da_);
  DMSetUp(da_);

  // Find out what are the corners of the subdomain
  DMDAGetCorners(da_, &firstX_, &firstY_, &firstZ_, &lengthX_, &lengthY_, &lengthZ_);

  // Current function to assign the limits
  createLimits(parameters, da_, limitsX_, limitsY_, limitsZ_);

  // Set the rank in the context. Necessary since PETSc declares the rank of the neighbors under
  // periodic conditions as the rank of the process itself. So the rank must be known to properly
  // set matrices and RHS vectors.
  int rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  // Set offsets to fix where the results of the pressure will be written in the flow field.
  if (firstX_ == 0) {
    offsetX_ = 1;
  } else {
    offsetX_ = 2;
  }

  if (firstY_ == 0) {
    offsetY_ = 1;
  } else {
    offsetY_ = 2;
  }

  if (firstZ_ == 0) {
    offsetZ_ = 1;
  } else {
    offsetZ_ = 2;
  }

  // Set a pointer to the limits in the context, so that they can be used by the function.
  ctx_.setLimits(limitsX_, limitsY_, limitsZ_);

  // Determine whether a process writes a boundary on the system.
  // Right now, it only depends on the position of the array. The identity of the neighbors will
  // become significant only in the iterators, where it will apply a boundary condition only if
  // the rank is invalid, so that there is no neighbor.
  ctx_.setAsBoundary = 0;
  if (parameters.parallel.indices[0] == 0) {
    ctx_.setAsBoundary += LEFT_WALL_BIT;
  }
  if (parameters.parallel.indices[0] == parameters.parallel.numProcessors[0] - 1) {
    ctx_.setAsBoundary += RIGHT_WALL_BIT;
  }
  if (parameters.parallel.indices[1] == 0) {
    ctx_.setAsBoundary += BOTTOM_WALL_BIT;
  }
  if (parameters.parallel.indices[1] == parameters.parallel.numProcessors[1] - 1) {
    ctx_.setAsBoundary += TOP_WALL_BIT;
  }
  if (parameters.parallel.indices[2] == 0) {
    ctx_.setAsBoundary += FRONT_WALL_BIT;
  }
  if (parameters.parallel.indices[2] == parameters.parallel.numProcessors[2] - 1) {
    ctx_.setAsBoundary += BACK_WALL_BIT;
  }

  // Set displacements to deal with periodic boundaries if necessary.
  // If the boundary is periodic, it will take information from positions beyond the ghost cells,
  // since they are used only for parallel communication. Otherwise, PETSc deals with
  // communication of the pressure.
  int Nx = parameters.geometry.sizeX + 2;
  int Ny = parameters.geometry.sizeY + 2;
  int Nz = parameters.geometry.sizeZ + 2;

  if (parameters.walls.typeLeft == PERIODIC) {
    ctx_.displacement[0] = -2;
    ctx_.displacement[1] = Nx + 1;
  } else {
    ctx_.displacement[0] = 1;
    ctx_.displacement[1] = Nx - 2;
  }

  if (parameters.walls.typeBottom == PERIODIC) {
    ctx_.displacement[2] = -2;
    ctx_.displacement[3] = Ny + 1;
  } else {
    ctx_.displacement[2] = 1;
    ctx_.displacement[3] = Ny - 2;
  }

  if (parameters.walls.typeFront == PERIODIC) {
    ctx_.displacement[4] = -2;
    ctx_.displacement[5] = Nz + 1;
  } else {
    ctx_.displacement[4] = 1;
    ctx_.displacement[5] = Nz - 2;
  }

  DMCreateGlobalVector(da_, &x_);
  KSPSetDM(ksp_, da_);
  KSPSetComputeOperators(ksp_, computeMatrix, &ctx_);

  KSPSetType(ksp_, KSPFGMRES);

  int commSize;
  MPI_Comm_size(PETSC_COMM_WORLD, &commSize);

  if (commSize == 1) {
    // If serial
    PCSetType(pc_, PCILU);
    PCFactorSetLevels(pc_, 1);
    KSPSetPC(ksp_, pc_);
  } else {
    // If parallel
    PCSetType(pc_, PCASM);
    KSPSetPC(ksp_, pc_);
  }

  KSPSetFromOptions(ksp_);
  KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE);
  KSPSetUp(ksp_);

  // From here we can change sub_ksp if necessary
  // that has to be done after setup. The other solvers above
  // can be changed before setup with KSPSetFromOptions.

  if (commSize > 1) {
    KSP* subksp;
    PC   subpc;

    PCASMGetSubKSP(pc_, NULL, NULL, &subksp);
    KSPGetPC(subksp[0], &subpc);

    PetscBool has_fl;
    PetscBool has_sub_type;
    PetscOptionsHasName(NULL, NULL, "-sub_pc_factor_levels", &has_fl);
    PetscOptionsHasName(NULL, NULL, "-sub_pc_type", &has_sub_type);
    if (!(has_sub_type))
      PCSetType(subpc, PCILU);
    if (!(has_fl))
      PCFactorSetLevels(subpc, 1);

    KSPSetUp(ksp_);
  }
}

void Solvers::PetscSolver::solve() {
  ScalarField& pressure = flowField_.getPressure();

  if (parameters_.geometry.dim == 2) {
    KSPSetComputeRHS(ksp_, computeRHS2D, &ctx_);
    KSPSetComputeOperators(ksp_, computeMatrix2D, &ctx_);
    KSPSolve(ksp_, PETSC_NULLPTR, x_);

    // Then extract the information
    PetscScalar** array;
    DMDAVecGetArray(da_, x_, &array);

    for (int j = firstY_; j < firstY_ + lengthY_; j++) {
      for (int i = firstX_; i < firstX_ + lengthX_; i++) {
        pressure.getScalar(i - firstX_ + offsetX_, j - firstY_ + offsetY_) = array[j][i];
      }
    }
    DMDAVecRestoreArray(da_, x_, &array);
  } else if (parameters_.geometry.dim == 3) {
    KSPSetComputeRHS(ksp_, computeRHS3D, &ctx_);
    KSPSetComputeOperators(ksp_, computeMatrix3D, &ctx_);
    KSPSolve(ksp_, PETSC_NULLPTR, x_);

    // Then extract the information
    PetscScalar*** array;
    DMDAVecGetArray(da_, x_, &array);

    for (int k = firstZ_; k < firstZ_ + lengthZ_; k++) {
      for (int j = firstY_; j < firstY_ + lengthY_; j++) {
        for (int i = firstX_; i < firstX_ + lengthX_; i++) {
          pressure.getScalar(i - firstX_ + offsetX_, j - firstY_ + offsetY_, k - firstZ_ + offsetZ_) = array[k][j][i];
        }
      }
    }
    DMDAVecRestoreArray(da_, x_, &array);
  }
}

PetscErrorCode computeMatrix2D([[maybe_unused]] KSP ksp, Mat A, [[maybe_unused]] Mat pc, void* ctx) {
  Solvers::PetscUserCtx* context    = static_cast<Solvers::PetscUserCtx*>(ctx);
  Parameters&            parameters = context->getParameters();

  IntScalarField& flags = context->getFlowField().getFlags();

  int *limitsX, *limitsY, *limitsZ;
  context->getLimits(&limitsX, &limitsY, &limitsZ);

  PetscScalar stencilValues[5];
  MatStencil  row, column[5];

  PetscInt i, j, Nx, Ny;

  Nx = parameters.geometry.sizeX + 2;
  Ny = parameters.geometry.sizeY + 2;

  spdlog::trace("Limits: {}, {}, {}, {}", limitsX[0], limitsX[1], limitsY[0], limitsY[1]);
  // Loop for inner nodes
  for (j = limitsY[0]; j < limitsY[1]; j++) {
    for (i = limitsX[0]; i < limitsX[1]; i++) {

      row.i = i;
      row.j = j;
      // Convert matrix to Cartesian grid coordinates
      const int cellIndexX = i - limitsX[0] + 2;
      const int cellIndexY = j - limitsY[0] + 2;

      const int obstacle = flags.getValue(cellIndexX, cellIndexY);

      if ((obstacle & OBSTACLE_SELF) == 0) { // If we have a fluid cell
        const RealType dx_0  = parameters.meshsize->getDx(cellIndexX, cellIndexY);
        const RealType dx_M1 = parameters.meshsize->getDx(cellIndexX - 1, cellIndexY);
        const RealType dx_P1 = parameters.meshsize->getDx(cellIndexX + 1, cellIndexY);
        const RealType dy_0  = parameters.meshsize->getDy(cellIndexX, cellIndexY);
        const RealType dy_M1 = parameters.meshsize->getDy(cellIndexX, cellIndexY - 1);
        const RealType dy_P1 = parameters.meshsize->getDy(cellIndexX, cellIndexY + 1);

        const RealType dx_L  = 0.5 * (dx_0 + dx_M1);
        const RealType dx_R  = 0.5 * (dx_0 + dx_P1);
        const RealType dx_Bo = 0.5 * (dy_0 + dy_M1);
        const RealType dx_T  = 0.5 * (dy_0 + dy_P1);

        // Definition of values: set general formulation for laplace operator here, based on arbitrary meshsizes.
        stencilValues[1] = 2.0 / (dx_L * (dx_L + dx_R));                // Left
        stencilValues[0] = 2.0 / (dx_R * (dx_L + dx_R));                // Right
        stencilValues[3] = 2.0 / (dx_T * (dx_T + dx_Bo));               // Top
        stencilValues[4] = 2.0 / (dx_Bo * (dx_T + dx_Bo));              // Bottom
        stencilValues[2] = -2.0 / (dx_R * dx_L) - 2.0 / (dx_T * dx_Bo); // Center

        // Definition of positions. Order must correspond to values.
        column[0].i = i + 1;
        column[0].j = j;
        column[1].i = i - 1;
        column[1].j = j;
        column[2].i = i;
        column[2].j = j;
        column[3].i = i;
        column[3].j = j + 1;
        column[4].i = i;
        column[4].j = j - 1;

        MatSetValuesStencil(A, 1, &row, 5, column, stencilValues, INSERT_VALUES);
      } else if (obstacle != OBSTACLE_SELF + OBSTACLE_LEFT + OBSTACLE_RIGHT + OBSTACLE_TOP + OBSTACLE_BOTTOM) { // Not
                                                                                                                // fluid,
                                                                                                                // but
                                                                                                                // fluid
                                                                                                                // somewhere
                                                                                                                // around.
        int counter      = 0; // This will contain how many neighbours are fluid
        int counterFluid = 0;
        // TODO: variable meshwidth might have to be considered
        if ((obstacle & OBSTACLE_LEFT) == 0) { // If there is fluid to the left
          stencilValues[counter] = 1.0;
          column[counter].i      = i - 1;
          column[counter].j      = j;
          counter++; // We have just identified a fuid cell and prepared to average
          counterFluid++;
        } else {
          stencilValues[counter] = 0.0;
          column[counter].i      = i - 1;
          column[counter].j      = j;
          counter++;
        }
        if ((obstacle & OBSTACLE_RIGHT) == 0) {
          stencilValues[counter] = 1.0;
          column[counter].i      = i + 1;
          column[counter].j      = j;
          counter++;
          counterFluid++;
        } else {
          stencilValues[counter] = 0.0;
          column[counter].i      = i + 1;
          column[counter].j      = j;
          counter++;
        }
        if ((obstacle & OBSTACLE_BOTTOM) == 0) {
          stencilValues[counter] = 1.0;
          column[counter].i      = i;
          column[counter].j      = j - 1;
          counter++;
          counterFluid++;
        } else {
          stencilValues[counter] = 0.0;
          column[counter].i      = i;
          column[counter].j      = j - 1;
          counter++;
        }
        if ((obstacle & OBSTACLE_TOP) == 0) {
          stencilValues[counter] = 1.0;
          column[counter].i      = i;
          column[counter].j      = j + 1;
          counter++;
          counterFluid++;
        } else {
          stencilValues[counter] = 0.0;
          column[counter].i      = i;
          column[counter].j      = j + 1;
          counter++;
        }

        // A column for the cell itself
        stencilValues[counter] = -static_cast<PetscScalar>(counterFluid);
        column[counter].i      = i;
        column[counter].j      = j;

        // Once we identified how many fluid cells are around and set columns for each, we
        // enter the row into the matrix.

        MatSetValuesStencil(A, 1, &row, 5, column, stencilValues, INSERT_VALUES);
      } else { // The remaining possibility is that the cell is obstacle surrounded
               // by more obstacle cells.
        // Here, we just add an equation to set the value according to the right hand side.
        stencilValues[1] = 0.0; // Left
        stencilValues[0] = 0.0; // Right
        stencilValues[3] = 0.0; // Top
        stencilValues[4] = 0.0; // Bottom
        stencilValues[2] = 1.0; // Center

        // Definition of positions. Order must correspond to values.
        column[0].i = i + 1;
        column[0].j = j;
        column[1].i = i - 1;
        column[1].j = j;
        column[2].i = i;
        column[2].j = j;
        column[3].i = i;
        column[3].j = j + 1;
        column[4].i = i;
        column[4].j = j - 1;
        MatSetValuesStencil(A, 1, &row, 5, column, stencilValues, INSERT_VALUES);
      }
    }
  }

  // Left wall
  if (context->setAsBoundary & LEFT_WALL_BIT) {
    for (j = limitsY[0]; j < limitsY[1]; j++) {
      column[0].i = 0;
      column[0].j = j;
      column[1].i = context->displacement[0];
      column[1].j = j;
      row.i       = 0;
      row.j       = j;
      if (parameters.walls.typeLeft == DIRICHLET) { // If Dirichlet velocity boundary conditions
        // Therefore, Neumann in the pressure
        stencilValues[0] = 1;
        stencilValues[1] = -1;
      } else if (parameters.walls.typeLeft == NEUMANN) { // Neumann velocity boundary conditions
        stencilValues[0] = 0.5;
        stencilValues[1] = 0.5;
      }
      MatSetValuesStencil(A, 1, &row, 2, column, stencilValues, INSERT_VALUES);
    }
  }

  // Right wall
  if (context->setAsBoundary & RIGHT_WALL_BIT) {
    for (j = limitsY[0]; j < limitsY[1]; j++) {
      column[0].i = Nx - 1;
      column[0].j = j;
      column[1].i = context->displacement[1];
      column[1].j = j;
      row.i       = Nx - 1;
      row.j       = j;
      if (parameters.walls.typeRight == DIRICHLET) {
        stencilValues[0] = 1;
        stencilValues[1] = -1;
      } else if (parameters.walls.typeRight == NEUMANN) {
        stencilValues[0] = 0.5;
        stencilValues[1] = 0.5;
      }
      MatSetValuesStencil(A, 1, &row, 2, column, stencilValues, INSERT_VALUES);
    }
  }

  // Bottom wall
  if (context->setAsBoundary & BOTTOM_WALL_BIT) {
    for (i = limitsX[0]; i < limitsX[1]; i++) {
      column[0].i = i;
      column[0].j = 0;
      column[1].i = i;
      column[1].j = context->displacement[2];
      row.i       = i;
      row.j       = 0;
      if (parameters.walls.typeBottom == DIRICHLET) {
        stencilValues[0] = 1;
        stencilValues[1] = -1;
      } else if (parameters.walls.typeBottom == NEUMANN) {
        stencilValues[0] = 0.5;
        stencilValues[1] = 0.5;
      }
      MatSetValuesStencil(A, 1, &row, 2, column, stencilValues, INSERT_VALUES);
    }
  }

  // Top wall
  if (context->setAsBoundary & TOP_WALL_BIT) {
    for (i = limitsX[0]; i < limitsX[1]; i++) {
      column[0].i = i;
      column[0].j = Ny - 1;
      column[1].i = i;
      column[1].j = context->displacement[3];
      row.i       = i;
      row.j       = Ny - 1;
      if (parameters.walls.typeTop == DIRICHLET) {
        stencilValues[0] = 1;
        stencilValues[1] = -1;
      } else if (parameters.walls.typeTop == NEUMANN) {
        stencilValues[0] = 0.5;
        stencilValues[1] = 0.5;
      }
      MatSetValuesStencil(A, 1, &row, 2, column, stencilValues, INSERT_VALUES);
    }
  }

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  MatNullSpace nullspace;
  MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace);
  MatSetNullSpace(A, nullspace);
  MatNullSpaceDestroy(&nullspace);

  return 0;
}

PetscErrorCode computeMatrix3D([[maybe_unused]] KSP ksp, Mat A, [[maybe_unused]] Mat pc, void* ctx) {
  Solvers::PetscUserCtx* context    = static_cast<Solvers::PetscUserCtx*>(ctx);
  Parameters&            parameters = context->getParameters();

  IntScalarField& flags = context->getFlowField().getFlags();

  int *limitsX, *limitsY, *limitsZ;
  context->getLimits(&limitsX, &limitsY, &limitsZ);

  PetscScalar stencilValues[7];
  MatStencil  row, column[7];

  PetscInt i, j, k, Nx, Ny, Nz;

  Nx = parameters.geometry.sizeX + 2;
  Ny = parameters.geometry.sizeY + 2;
  Nz = parameters.geometry.sizeZ + 2;

  // Loop for inner nodes
  spdlog::trace(
    "Limits: {}, {}, {}, {}, {}, {}", limitsX[0], limitsX[1], limitsY[0], limitsY[1], limitsZ[0], limitsZ[1]
  );
  for (k = limitsZ[0]; k < limitsZ[1]; k++) {
    for (j = limitsY[0]; j < limitsY[1]; j++) {
      for (i = limitsX[0]; i < limitsX[1]; i++) {
        row.i = i;
        row.j = j, row.k = k;

        const int cellIndexX = i - limitsX[0] + 2;
        const int cellIndexY = j - limitsY[0] + 2;
        const int cellIndexZ = k - limitsZ[0] + 2;
        const int obstacle   = flags.getValue(cellIndexX, cellIndexY, cellIndexZ);

        if ((obstacle & OBSTACLE_SELF) == 0) { // If the cell is fluid
          const RealType dx_0  = parameters.meshsize->getDx(cellIndexX, cellIndexY, cellIndexZ);
          const RealType dx_M1 = parameters.meshsize->getDx(cellIndexX - 1, cellIndexY, cellIndexZ);
          const RealType dx_P1 = parameters.meshsize->getDx(cellIndexX + 1, cellIndexY, cellIndexZ);
          const RealType dy_0  = parameters.meshsize->getDy(cellIndexX, cellIndexY, cellIndexZ);
          const RealType dy_M1 = parameters.meshsize->getDy(cellIndexX, cellIndexY - 1, cellIndexZ);
          const RealType dy_P1 = parameters.meshsize->getDy(cellIndexX, cellIndexY + 1, cellIndexZ);
          const RealType dz_0  = parameters.meshsize->getDz(cellIndexX, cellIndexY, cellIndexZ);
          const RealType dz_M1 = parameters.meshsize->getDz(cellIndexX, cellIndexY, cellIndexZ - 1);
          const RealType dz_P1 = parameters.meshsize->getDz(cellIndexX, cellIndexY, cellIndexZ + 1);

          const RealType dx_L  = 0.5 * (dx_0 + dx_M1);
          const RealType dx_R  = 0.5 * (dx_0 + dx_P1);
          const RealType dx_Bo = 0.5 * (dy_0 + dy_M1);
          const RealType dx_T  = 0.5 * (dy_0 + dy_P1);
          const RealType dx_F  = 0.5 * (dz_0 + dz_M1);
          const RealType dx_B  = 0.5 * (dz_0 + dz_P1);

          // Definition of values
          stencilValues[1] = 2.0 / (dx_L * (dx_L + dx_R));                                      // Left
          stencilValues[0] = 2.0 / (dx_R * (dx_L + dx_R));                                      // Right
          stencilValues[2] = 2.0 / (dx_T * (dx_T + dx_Bo));                                     // Top
          stencilValues[3] = 2.0 / (dx_Bo * (dx_T + dx_Bo));                                    // Bottom
          stencilValues[4] = 2.0 / (dx_B * (dx_B + dx_F));                                      // Back
          stencilValues[5] = 2.0 / (dx_F * (dx_B + dx_F));                                      // Front
          stencilValues[6] = -2.0 / (dx_R * dx_L) - 2.0 / (dx_T * dx_Bo) - 2.0 / (dx_F * dx_B); // Center

          // Definition of positions. Order must correspond to values.
          column[0].i = i + 1;
          column[0].j = j;
          column[0].k = k;
          column[1].i = i - 1;
          column[1].j = j;
          column[1].k = k;
          column[2].i = i;
          column[2].j = j + 1;
          column[2].k = k;
          column[3].i = i;
          column[3].j = j - 1;
          column[3].k = k;
          column[4].i = i;
          column[4].j = j;
          column[4].k = k + 1;
          column[5].i = i;
          column[5].j = j;
          column[5].k = k - 1;
          column[6].i = i;
          column[6].j = j;
          column[6].k = k;

          MatSetValuesStencil(A, 1, &row, 7, column, stencilValues, INSERT_VALUES);
        } else if (obstacle != 127) { // If non-fluid and still not completely surounded
          int counter      = 0;
          int counterFluid = 0;
          if ((obstacle & OBSTACLE_LEFT) == 0) { // If there's fluid to the left
            stencilValues[counter] = 1;
            column[counter].i      = i - 1;
            column[counter].j = j, column[counter].k = k;
            counter++;
            counterFluid++;
          } else {
            stencilValues[counter] = 0.0;
            column[counter].i      = i - 1;
            column[counter].j = j, column[counter].k = k;
            counter++;
          }
          if ((obstacle & OBSTACLE_RIGHT) == 0) {
            stencilValues[counter] = 1;
            column[counter].i      = i + 1;
            column[counter].j = j, column[counter].k = k;
            counter++;
            counterFluid++;
          } else {
            stencilValues[counter] = 0.0;
            column[counter].i      = i + 1;
            column[counter].j = j, column[counter].k = k;
            counter++;
          }
          if ((obstacle & OBSTACLE_BOTTOM) == 0) {
            stencilValues[counter] = 1;
            column[counter].i      = i;
            column[counter].j = j - 1, column[counter].k = k;
            counter++;
            counterFluid++;
          } else {
            stencilValues[counter] = 0.0;
            column[counter].i      = i;
            column[counter].j = j - 1, column[counter].k = k;
            counter++;
          }
          if ((obstacle & OBSTACLE_TOP) == 0) {
            stencilValues[counter] = 1;
            column[counter].i      = i;
            column[counter].j = j + 1, column[counter].k = k;
            counter++;
            counterFluid++;
          } else {
            stencilValues[counter] = 0.0;
            column[counter].i      = i;
            column[counter].j = j + 1, column[counter].k = k;
            counter++;
          }
          if ((obstacle & OBSTACLE_FRONT) == 0) {
            stencilValues[counter] = 1;
            column[counter].i      = i;
            column[counter].j = j, column[counter].k = k - 1;
            counter++;
            counterFluid++;
          } else {
            stencilValues[counter] = 0.0;
            column[counter].i      = i;
            column[counter].j = j, column[counter].k = k - 1;
            counter++;
          }
          if ((obstacle & OBSTACLE_BACK) == 0) {
            stencilValues[counter] = 1;
            column[counter].i      = i;
            column[counter].j = j, column[counter].k = k + 1;
            counter++;
            counterFluid++;
          } else {
            stencilValues[counter] = 0.0;
            column[counter].i      = i;
            column[counter].j = j, column[counter].k = k + 1;
            counter++;
          }

          // Now set a line for the element itself
          stencilValues[counter] = static_cast<PetscScalar>(-counterFluid);
          column[counter].i      = i;
          column[counter].j = j, column[counter].k = k;

          MatSetValuesStencil(A, 1, &row, counter + 1, column, stencilValues, INSERT_VALUES);
        } else { // The remaining possibility is that the cell is obstacle surrounded
                 // by more obstacle cells.
          // Here, we just add an equation to set the value according to the right hand side.
          stencilValues[1] = 0.0; // Left
          stencilValues[0] = 0.0; // Right
          stencilValues[3] = 0.0; // Top
          stencilValues[4] = 0.0; // Bottom
          stencilValues[5] = 0.0; // Front
          stencilValues[5] = 0.0; // Back
          stencilValues[2] = 1.0; // Center

          // Definition of positions. Order must correspond to values
          column[0].i = i + 1;
          column[0].j = j;
          column[0].k = k;
          column[1].i = i - 1;
          column[1].j = j;
          column[1].k = k;
          column[2].i = i;
          column[2].j = j;
          column[2].k = k;
          column[3].i = i;
          column[3].j = j + 1;
          column[3].k = k;
          column[4].i = i;
          column[4].j = j - 1;
          column[4].k = k;
          column[5].i = i;
          column[5].j = j;
          column[5].k = k + 1;
          column[6].i = i;
          column[6].j = j;
          column[6].k = k - 1;
          MatSetValuesStencil(A, 1, &row, 7, column, stencilValues, INSERT_VALUES);
        }
      }
    }
  }

  // Left wall
  if (context->setAsBoundary & LEFT_WALL_BIT) {
    for (j = limitsY[0]; j < limitsY[1]; j++) {
      for (k = limitsZ[0]; k < limitsZ[1]; k++) {
        column[0].i = 0;
        column[0].j = j;
        column[0].k = k;
        column[1].i = context->displacement[0];
        column[1].j = j;
        column[1].k = k;
        row.i       = 0;
        row.j       = j;
        row.k       = k;
        if (parameters.walls.typeLeft == DIRICHLET) {
          stencilValues[0] = 1;
          stencilValues[1] = -1;
        } else if (parameters.walls.typeLeft == NEUMANN) {
          stencilValues[0] = 0.5;
          stencilValues[1] = 0.5;
        }
        MatSetValuesStencil(A, 1, &row, 2, column, stencilValues, INSERT_VALUES);
      }
    }
  }

  // Right wall
  if (context->setAsBoundary & RIGHT_WALL_BIT) {
    for (j = limitsY[0]; j < limitsY[1]; j++) {
      for (k = limitsZ[0]; k < limitsZ[1]; k++) {
        column[0].i = Nx - 1;
        column[0].j = j;
        column[0].k = k;
        column[1].i = context->displacement[1];
        column[1].j = j;
        column[1].k = k;
        row.i       = Nx - 1;
        row.j       = j;
        row.k       = k;
        if (parameters.walls.typeRight == DIRICHLET) {
          stencilValues[0] = 1;
          stencilValues[1] = -1;
        } else if (parameters.walls.typeRight == NEUMANN) {
          stencilValues[0] = 0.5;
          stencilValues[1] = 0.5;
        }
        MatSetValuesStencil(A, 1, &row, 2, column, stencilValues, INSERT_VALUES);
      }
    }
  }

  // Bottom wall
  if (context->setAsBoundary & BOTTOM_WALL_BIT) {
    for (i = limitsX[0]; i < limitsX[1]; i++) {
      for (k = limitsZ[0]; k < limitsZ[1]; k++) {
        column[0].i = i;
        column[0].j = 0;
        column[0].k = k;
        column[1].i = i;
        column[1].j = context->displacement[2];
        column[1].k = k;
        row.i       = i;
        row.j       = 0;
        row.k       = k;
        if (parameters.walls.typeBottom == DIRICHLET) {
          stencilValues[0] = 1;
          stencilValues[1] = -1;
        } else if (parameters.walls.typeBottom == NEUMANN) {
          stencilValues[0] = 0.5;
          stencilValues[1] = 0.5;
        }
        MatSetValuesStencil(A, 1, &row, 2, column, stencilValues, INSERT_VALUES);
      }
    }
  }

  // Top wall
  if (context->setAsBoundary & TOP_WALL_BIT) {
    for (i = limitsX[0]; i < limitsX[1]; i++) {
      for (k = limitsZ[0]; k < limitsZ[1]; k++) {
        column[0].i = i;
        column[0].j = Ny - 1;
        column[0].k = k;
        column[1].i = i;
        column[1].j = context->displacement[3];
        column[1].k = k;
        row.i       = i;
        row.j       = Ny - 1;
        row.k       = k;
        if (parameters.walls.typeTop == DIRICHLET) {
          stencilValues[0] = 1;
          stencilValues[1] = -1;
        } else if (parameters.walls.typeTop == NEUMANN) {
          stencilValues[0] = 0.5;
          stencilValues[1] = 0.5;
        }
        MatSetValuesStencil(A, 1, &row, 2, column, stencilValues, INSERT_VALUES);
      }
    }
  }

  // Front wall
  if (context->setAsBoundary & FRONT_WALL_BIT) {
    for (i = limitsX[0]; i < limitsX[1]; i++) {
      for (j = limitsY[0]; j < limitsY[1]; j++) {
        column[0].i = i;
        column[0].j = j;
        column[0].k = 0;
        column[1].i = i;
        column[1].j = j;
        column[1].k = context->displacement[4];
        row.i       = i;
        row.j       = j;
        row.k       = 0;
        if (parameters.walls.typeFront == DIRICHLET) {
          stencilValues[0] = 1;
          stencilValues[1] = -1;
        } else if (parameters.walls.typeFront == NEUMANN) {
          stencilValues[0] = 0.5;
          stencilValues[1] = 0.5;
        }
        MatSetValuesStencil(A, 1, &row, 2, column, stencilValues, INSERT_VALUES);
      }
    }
  }

  // Back wall
  if (context->setAsBoundary & BACK_WALL_BIT) {
    for (i = limitsX[0]; i < limitsX[1]; i++) {
      for (j = limitsY[0]; j < limitsY[1]; j++) {
        column[0].i = i;
        column[0].j = j;
        column[0].k = Nz - 1;
        column[1].i = i;
        column[1].j = j;
        column[1].k = context->displacement[5];
        row.i       = i;
        row.j       = j;
        row.k       = Nz - 1;
        if (parameters.walls.typeBack == DIRICHLET) {
          stencilValues[0] = 1;
          stencilValues[1] = -1;
        } else if (parameters.walls.typeBack == NEUMANN) {
          stencilValues[0] = 0.5;
          stencilValues[1] = 0.5;
        }
        MatSetValuesStencil(A, 1, &row, 2, column, stencilValues, INSERT_VALUES);
      }
    }
  }

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  MatNullSpace nullspace;
  MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace);
  MatSetNullSpace(A, nullspace);
  MatNullSpaceDestroy(&nullspace);

  return 0;
}

PetscErrorCode computeRHS2D(KSP ksp, Vec b, void* ctx) {
  FlowField&             flowField  = static_cast<Solvers::PetscUserCtx*>(ctx)->getFlowField();
  Parameters&            parameters = static_cast<Solvers::PetscUserCtx*>(ctx)->getParameters();
  Solvers::PetscUserCtx* context    = static_cast<Solvers::PetscUserCtx*>(ctx);

  int *limitsX, *limitsY, *limitsZ;
  static_cast<Solvers::PetscUserCtx*>(ctx)->getLimits(&limitsX, &limitsY, &limitsZ);

  IntScalarField& flags = context->getFlowField().getFlags();

  ScalarField& RHS = flowField.getRHS();

  PetscInt      i, j;
  PetscInt      Nx = parameters.geometry.sizeX + 2, Ny = parameters.geometry.sizeY + 2;
  PetscScalar** array;

  DM da;
  KSPGetDM(ksp, &da);

  DMDAVecGetArray(da, b, &array);

  // Iteration domains are going to be set and the values on the global boundary set when necessary
  // Check left wall
  if (context->setAsBoundary & LEFT_WALL_BIT) {
    if (parameters.simulation.scenario == "pressure-channel") {
      for (j = limitsY[0]; j < limitsY[1]; j++) {
        array[j][0] = RHS.getScalar(0, j);
      }
    } else {
      for (j = limitsY[0]; j < limitsY[1]; j++) {
        array[j][0] = 0;
      }
    }
  }

  // Check right wall
  if (context->setAsBoundary & RIGHT_WALL_BIT) {
    for (j = limitsY[0]; j < limitsY[1]; j++) {
      array[j][Nx - 1] = 0;
    }
  }

  // Check bottom wall
  if (context->setAsBoundary & BOTTOM_WALL_BIT) {
    for (i = limitsX[0]; i < limitsX[1]; i++) {
      array[0][i] = 0;
    }
  }

  // Check top wall
  if (context->setAsBoundary & TOP_WALL_BIT) {
    for (i = limitsX[0]; i < limitsX[1]; i++) {
      array[Ny - 1][i] = 0;
    }
  }

  // Fill the internal nodes. We already have the values.
  for (j = limitsY[0]; j < limitsY[1]; j++) {
    for (i = limitsX[0]; i < limitsX[1]; i++) {
      const int obstacle = flags.getValue(i - limitsX[0] + 2, j - limitsY[0] + 2);
      if ((obstacle & OBSTACLE_SELF) == 0) { // If this is a fluid cell
        array[j][i] = RHS.getScalar(i - limitsX[0] + 2, j - limitsY[0] + 2);
      } else {
        array[j][i] = 0.0;
      }
    }
  }

  DMDAVecRestoreArray(da, b, &array);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  return 0;
}

PetscErrorCode computeRHS3D(KSP ksp, Vec b, void* ctx) {
  FlowField&             flowField  = static_cast<Solvers::PetscUserCtx*>(ctx)->getFlowField();
  Parameters&            parameters = static_cast<Solvers::PetscUserCtx*>(ctx)->getParameters();
  ScalarField&           RHS        = flowField.getRHS();
  Solvers::PetscUserCtx* context    = static_cast<Solvers::PetscUserCtx*>(ctx);

  IntScalarField& flags = flowField.getFlags();

  int *limitsX, *limitsY, *limitsZ;
  static_cast<Solvers::PetscUserCtx*>(ctx)->getLimits(&limitsX, &limitsY, &limitsZ);

  PetscInt i, j, k;
  PetscInt Nx = parameters.geometry.sizeX + 2, Ny = parameters.geometry.sizeY + 2, Nz = parameters.geometry.sizeZ + 2;
  PetscScalar*** array;

  DM da;
  KSPGetDM(ksp, &da);

  DMDAVecGetArray(da, b, &array);

  // Notice that we're covering the whole surface, including corners and edges
  // Also, the actual value is taking from the parameters.

  // Left wall
  if (context->setAsBoundary & LEFT_WALL_BIT) {
    if (parameters.simulation.scenario == "pressure-channel") {
      for (k = limitsZ[0]; k < limitsZ[1]; k++) {
        for (j = limitsY[0]; j < limitsY[1]; j++) {
          array[k][j][0] = RHS.getScalar(0, j, k);
        }
      }
    } else {
      for (k = limitsZ[0]; k < limitsZ[1]; k++) {
        for (j = limitsY[0]; j < limitsY[1]; j++) {
          array[k][j][0] = 0;
        }
      }
    }
  }

  // Right wall
  if (context->setAsBoundary & RIGHT_WALL_BIT) {
    for (k = limitsZ[0]; k < limitsZ[1]; k++) {
      for (j = limitsY[0]; j < limitsY[1]; j++) {
        array[k][j][Nx - 1] = 0;
      }
    }
  }

  // Bottom wall
  if (context->setAsBoundary & BOTTOM_WALL_BIT) {
    for (k = limitsZ[0]; k < limitsZ[1]; k++) {
      for (i = limitsX[0]; i < limitsX[1]; i++) {
        array[k][0][i] = 0;
      }
    }
  }

  // Top wall
  if (context->setAsBoundary & TOP_WALL_BIT) {
    for (k = limitsZ[0]; k < limitsZ[1]; k++) {
      for (i = limitsX[0]; i < limitsX[1]; i++) {
        array[k][Ny - 1][i] = 0;
      }
    }
  }

  // Front wall
  if (context->setAsBoundary & FRONT_WALL_BIT) {
    for (j = limitsY[0]; j < limitsY[1]; j++) {
      for (i = limitsX[0]; i < limitsX[1]; i++) {
        array[0][j][i] = 0;
      }
    }
  }

  // Back wall
  if (context->setAsBoundary & BACK_WALL_BIT) {
    for (j = limitsY[0]; j < limitsY[1]; j++) {
      for (i = limitsX[0]; i < limitsX[1]; i++) {
        array[Nz - 1][j][i] = 0;
      }
    }
  }

  // Fill the internal nodes. We already have the values.
  for (k = limitsZ[0]; k < limitsZ[1]; k++) {
    for (j = limitsY[0]; j < limitsY[1]; j++) {
      for (i = limitsX[0]; i < limitsX[1]; i++) {
        const int obstacle = flags.getValue(i - limitsX[0] + 2, j - limitsY[0] + 2, k - limitsZ[0] + 2);
        if ((obstacle & OBSTACLE_SELF) == 0) { // If this is a fluid cell
          array[k][j][i] = RHS.getScalar(i - limitsX[0] + 2, j - limitsY[0] + 2, k - limitsZ[0] + 2);
        } else {
          array[k][j][i] = 0.0;
        }
      }
    }
  }

  DMDAVecRestoreArray(da, b, &array);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);

  return 0;
}

const DM& Solvers::PetscSolver::getGrid() const { return da_; }

void Solvers::PetscUserCtx::setLimits(int* limitsX, int* limitsY, int* limitsZ) {
  limitsX_ = limitsX;
  limitsY_ = limitsY;
  limitsZ_ = limitsZ;
}

void Solvers::PetscUserCtx::getLimits(int** limitsX, int** limitsY, int** limitsZ) {
  *limitsX = limitsX_;
  *limitsY = limitsY_;
  *limitsZ = limitsZ_;
}

void Solvers::PetscUserCtx::setRank(int rank) { rank_ = rank; }

int Solvers::PetscUserCtx::getRank() const { return rank_; }

void Solvers::PetscSolver::reInitMatrix() {
  spdlog::info("Reinit the matrix");
  if (parameters_.geometry.dim == 2) {
    KSPSetComputeOperators(ksp_, computeMatrix2D, &ctx_);
  } else {
    KSPSetComputeOperators(ksp_, computeMatrix3D, &ctx_);
  }
}

#endif
