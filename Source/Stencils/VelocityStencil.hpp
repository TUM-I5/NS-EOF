#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /** Stencil to compute the velocity once the pressure has been found.
   */
  class VelocityStencil: public FieldStencil<FlowField> {
  public:
    VelocityStencil(const Parameters& parameters);
    ~VelocityStencil() override = default;

    void apply(FlowField& flowField, int i, int j) override;
    void apply(FlowField& flowField, int i, int j, int k) override;
  };

} // namespace Stencils
