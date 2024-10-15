#pragma once

#include "BoundaryStencil.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /** Computes the maximum value of max(velocity)/meshsize for all grid cells.
   *
   * Originally, one would compute the max. velocity only and adapt it with the meshsize afterwards.
   * This, however, becomes inconsistent when dealing with non-uniform, e.g. stretched, meshes, since
   * the meshsize may be different for every grid cell. We therefore determine the max(velocity)/meshsize
   * and synchronise this value over whole computational domain.
   */
  class MaxUStencil: public FieldStencil<FlowField>, public BoundaryStencil<FlowField> {
  private:
    RealType maxValues_[3]; //! Stores the maximum module of every component

    /** Sets the maximum value arrays to the value of the cell if it surpasses the current one.
     *
     * 2D version of the function
     * @param flowField Flow field
     * @param i Position in the X direction.
     * @param j Position in the Y direction.
     */
    void cellMaxValue(FlowField& flowField, int i, int j);

    /** Sets the maximum value arrays to the value of the cell if it surpasses the current one.
     *
     * 3D version of the function
     * @param flowField Flow field
     * @param i Position in the X direction.
     * @param j Position in the Y direction.
     * @param k Position in the Z direction.
     */
    void cellMaxValue(FlowField& flowField, int i, int j, int k);

  public:
    MaxUStencil(const Parameters& parameters);
    ~MaxUStencil() override = default;

    void apply(FlowField& flowField, int i, int j) override;
    void apply(FlowField& flowField, int i, int j, int k) override;

    void applyLeftWall(FlowField& flowField, int i, int j) override;
    void applyRightWall(FlowField& flowField, int i, int j) override;
    void applyBottomWall(FlowField& flowField, int i, int j) override;
    void applyTopWall(FlowField& flowField, int i, int j) override;

    void applyLeftWall(FlowField& flowField, int i, int j, int k) override;
    void applyRightWall(FlowField& flowField, int i, int j, int k) override;
    void applyBottomWall(FlowField& flowField, int i, int j, int k) override;
    void applyTopWall(FlowField& flowField, int i, int j, int k) override;
    void applyFrontWall(FlowField& flowField, int i, int j, int k) override;
    void applyBackWall(FlowField& flowField, int i, int j, int k) override;

    /** Resets the maximum values to zero before computing the timestep.
     */
    void reset();

    /** Returns the array with the maximum modules of the components of the velocity,
     *  divided by the respective local meshsize.
     */
    const RealType* getMaxValues() const;
  };

} // namespace Stencils
