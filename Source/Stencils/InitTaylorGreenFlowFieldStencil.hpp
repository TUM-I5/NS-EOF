#pragma once

#include "Definitions.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /** Stencil for initialising the taylor-green vortex flow in a periodic domain.
   */
  class InitTaylorGreenFlowFieldStencil: public FieldStencil<FlowField> {
  private:
    const RealType        pi2_;
    const RealType* const domainSize_;

    RealType* initializeDomainSize(const Parameters& parameters) const;

    /** from the local grid coordinates i, j, k, computes the global coordinates of the current cell and initialises
     *  the velocity field correspondingly.
     */
    void computeGlobalCoordinates(RealType* coords, int i, int j, int k = 0) const;

  public:
    InitTaylorGreenFlowFieldStencil(const Parameters& parameters);
    ~InitTaylorGreenFlowFieldStencil() override;

    void apply(FlowField& flowField, int i, int j) override;
    void apply(FlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils
