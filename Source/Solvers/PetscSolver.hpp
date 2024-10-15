#pragma once

#include "DataStructures.hpp"
#include "LinearSolver.hpp"

#ifdef ENABLE_PETSC

namespace Solvers {

  /**
   * A class to encapsulate information the Petsc builder functions.
   *  Petsc used so called context objects to give information to its routines.
   *  We need them to pass the flow field and the parameters in a single argument.
   */
  class PetscUserCtx {
  private:
    Parameters& parameters_;
    FlowField&  flowField_;

    int *limitsX_, *limitsY_, *limitsZ_;

    int rank_;

  public:
    PetscUserCtx(Parameters& parameters, FlowField& flowField);
    ~PetscUserCtx() = default;

    Parameters& getParameters();
    FlowField&  getFlowField();

    void setLimits(int* limitsX, int* limitsY, int* limitsZ);
    void getLimits(int** limitsX, int** limitsY, int** limitsZ);

    void setRank(int rank);
    int  getRank() const;

    unsigned char setAsBoundary;   // If set as boundary in the linear system. Use bits.
    int           displacement[6]; // Displacements for the boundary treatment
  };

  class PetscSolver: public LinearSolver {
  private:
    Vec x_;   //! Petsc vectors for solution and RHS
    DM  da_;  //! Topology manager
    KSP ksp_; //! Solver context
    PC  pc_;  //! Preconditioner

    PetscUserCtx ctx_; //! Capsule for Petsc builders

    // Indices for filling the matrices and right hand side
    int limitsX_[2], limitsY_[2], limitsZ_[2];

    PetscInt firstX_, lengthX_, firstY_, lengthY_, firstZ_, lengthZ_;

    // Additional variables used to determine where to write back the results
    int offsetX_, offsetY_, offsetZ_;

  public:
    PetscSolver(FlowField& flowField, Parameters& parameters);
    ~PetscSolver() override = default;

    void solve() override;

    // Reinit the matrix so that it uses the right flag field
    void reInitMatrix() override;

    const DM& getGrid() const;
  };

} // namespace Solvers

#endif
