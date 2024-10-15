#pragma once

#include "Parameters.hpp"

namespace Stencils {

  /** Stencil class
   *
   * Abstract class for the definition of stencils and operations on the grids
   */
  template <class FlowFieldType>
  class FieldStencil {
  protected:
    const Parameters& parameters_; //! Reference to the parameters

  public:
    FieldStencil(const Parameters& parameters):
      parameters_(parameters) {}

    virtual ~FieldStencil() = default;

    /** Performs the operation in 2D in a given position
     * @param flowField Flow field data
     * @param parameters Parameters of the problem
     * @param i Position in the x direction
     * @param j Position in the y direction
     */
    virtual void apply(FlowFieldType& flowField, int i, int j) = 0;

    /** Performs the operation in 3D in a given position
     * @param flowField Flow field data
     * @param parameters Parameters of the problem
     * @param i Position in the x direction
     * @param j Position in the y direction
     * @param k Position in the z direction
     */
    virtual void apply(FlowFieldType& flowField, int i, int j, int k) = 0;
  };

} // namespace Stencils
