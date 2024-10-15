#pragma once

#include "Parameters.hpp"

namespace Stencils {

  /** Interface for operations on the global (or parallel) boundary
   */
  template <class FlowFieldType>
  class BoundaryStencil {
  protected:
    const Parameters& parameters_;

  public:
    BoundaryStencil(const Parameters& parameters):
      parameters_(parameters) {}

    virtual ~BoundaryStencil() = default;

    /** Represents an operation in the left wall of a 2D domain.
     *
     * @param flowField State of the flow field
     * @param i Index in the x direction
     * @param j Index in the y direction
     */
    virtual void applyLeftWall(FlowFieldType& flowField, int i, int j) = 0;

    /** Represents an operation in the right wall of a 2D domain.
     *
     * @param flowField State of the flow field
     * @param i Index in the x direction
     * @param j Index in the y direction
     */
    virtual void applyRightWall(FlowFieldType& flowField, int i, int j) = 0;

    /** Represents an operation in the bottom wall of a 2D domain.
     *
     * @param flowField State of the flow field
     * @param i Index in the x direction
     * @param j Index in the y direction
     */
    virtual void applyBottomWall(FlowFieldType& flowField, int i, int j) = 0;

    /** Represents an operation in the top wall of a 2D domain.
     *
     * @param flowField State of the flow field
     * @param i Index in the x direction
     * @param j Index in the y direction
     */
    virtual void applyTopWall(FlowFieldType& flowField, int i, int j) = 0;

    /** Represents an operation in the left wall of a 3D domain.
     *
     * @param flowField State of the flow field
     * @param i Index in the x direction
     * @param j Index in the y direction
     * @param k Index in the z direction
     */
    virtual void applyLeftWall(FlowFieldType& flowField, int i, int j, int k) = 0;

    /** Represents an operation in the right wall of a 3D domain.
     *
     * @param flowField State of the flow field
     * @param i Index in the x direction
     * @param j Index in the y direction
     * @param k Index in the z direction
     */
    virtual void applyRightWall(FlowFieldType& flowField, int i, int j, int k) = 0;

    /** Represents an operation in the bottom wall of a 3D domain.
     *
     * @param flowField State of the flow field
     * @param i Index in the x direction
     * @param j Index in the y direction
     * @param k Index in the z direction
     */
    virtual void applyBottomWall(FlowFieldType& flowField, int i, int j, int k) = 0;

    /** Represents an operation in the top wall of a 3D domain.
     *
     * @param flowField State of the flow field
     * @param i Index in the x direction
     * @param j Index in the y direction
     * @param k Index in the z direction
     */
    virtual void applyTopWall(FlowFieldType& flowField, int i, int j, int k) = 0;

    /** Represents an operation in the front wall of a 3D domain.
     *
     * @param flowField State of the flow field
     * @param i Index in the x direction
     * @param j Index in the y direction
     * @param k Index in the z direction
     */
    virtual void applyFrontWall(FlowFieldType& flowField, int i, int j, int k) = 0;

    /** Represents an operation in the back wall of a 3D domain.
     *
     * @param flowField State of the flow field
     * @param i Index in the x direction
     * @param j Index in the y direction
     * @param k Index in the z direction
     */
    virtual void applyBackWall(FlowFieldType& flowField, int i, int j, int k) = 0;
  };

} // namespace Stencils
