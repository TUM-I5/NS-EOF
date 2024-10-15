#pragma once

#include "Parameters.hpp"

#include "Stencils/BoundaryStencil.hpp"
#include "Stencils/FieldStencil.hpp"

/** Iterator class
 *
 * Applies operations to a flow field
 */
template <class FlowFieldType>
class Iterator {
protected:
  FlowFieldType&    flowField_;
  const Parameters& parameters_;

public:
  Iterator(FlowFieldType& flowfield, const Parameters& parameters):
    flowField_(flowfield),
    parameters_(parameters) {}

  virtual ~Iterator() = default;

  /** Perform the stencil operation on inner, non-ghost cells.
   */
  virtual void iterate() = 0;
};

template <class FlowFieldType>
class FieldIterator: public Iterator<FlowFieldType> {
private:
  Stencils::FieldStencil<FlowFieldType>& stencil_;

  //@brief Define the iteration domain to include more or less layers
  // Added since the ability to select the iteration domain provides more flexibility
  //@{
  const int lowOffset_;
  const int highOffset_;
  //@}

public:
  FieldIterator(
    FlowFieldType&                         flowField,
    const Parameters&                      parameters,
    Stencils::FieldStencil<FlowFieldType>& stencil,
    int                                    lowOffset  = 0,
    int                                    highOffset = 0
  );

  virtual ~FieldIterator() override = default;

  /** Volume iteration over the field.
   *
   * Volume iteration. The stencil will be applied to all cells in the domain plus the upper
   * boundaries. Lower boundaries are not included.
   */
  virtual void iterate() override;
};

template <class FlowFieldType>
class GlobalBoundaryIterator: public Iterator<FlowFieldType> {
private:
  const int lowOffset_;
  const int highOffset_;

  // This iterator has a reference to a stencil for each side, and will call its methods.
  Stencils::BoundaryStencil<FlowFieldType>& leftWallStencil_;
  Stencils::BoundaryStencil<FlowFieldType>& rightWallStencil_;
  Stencils::BoundaryStencil<FlowFieldType>& bottomWallStencil_;
  Stencils::BoundaryStencil<FlowFieldType>& topWallStencil_;
  Stencils::BoundaryStencil<FlowFieldType>& frontWallStencil_;
  Stencils::BoundaryStencil<FlowFieldType>& backWallStencil_;

public:
  GlobalBoundaryIterator(
    FlowFieldType&                            flowField,
    const Parameters&                         parameters,
    Stencils::BoundaryStencil<FlowFieldType>& stencil,
    int                                       lowOffset  = 0,
    int                                       highOffset = 0
  );

  GlobalBoundaryIterator(
    FlowFieldType&                            flowField,
    const Parameters&                         parameters,
    Stencils::BoundaryStencil<FlowFieldType>& leftWallStencil,
    Stencils::BoundaryStencil<FlowFieldType>& rightWallStencil,
    Stencils::BoundaryStencil<FlowFieldType>& bottomWallStencil,
    Stencils::BoundaryStencil<FlowFieldType>& topWallStencil,
    int                                       lowOffset  = 0,
    int                                       highOffset = 0
  );

  GlobalBoundaryIterator(
    FlowFieldType&                            flowField,
    const Parameters&                         parameters,
    Stencils::BoundaryStencil<FlowFieldType>& leftWallStencil,
    Stencils::BoundaryStencil<FlowFieldType>& rightWallStencil,
    Stencils::BoundaryStencil<FlowFieldType>& bottomWallStencil,
    Stencils::BoundaryStencil<FlowFieldType>& topWallStencil,
    Stencils::BoundaryStencil<FlowFieldType>& frontWallStencil,
    Stencils::BoundaryStencil<FlowFieldType>& backWallStencil,
    int                                       lowOffset  = 0,
    int                                       highOffset = 0
  );

  virtual ~GlobalBoundaryIterator() override = default;

  /** Surface iterator
   *
   * Iterates on the boundary cells. Only upper corners and edges are iterated.
   */
  virtual void iterate() override;
};

template <class FlowFieldType>
class ParallelBoundaryIterator: public Iterator<FlowFieldType> {
private:
  Stencils::BoundaryStencil<FlowFieldType>& stencil_;

  const int lowOffset_;
  const int highOffset_;

public:
  ParallelBoundaryIterator(
    FlowFieldType&                            flowField,
    const Parameters&                         parameters,
    Stencils::BoundaryStencil<FlowFieldType>& stencil,
    int                                       lowOffset  = 0,
    int                                       highOffset = 0
  );

  virtual ~ParallelBoundaryIterator() override = default;

  virtual void iterate() override;
};

#include "Iterators.cpph"
