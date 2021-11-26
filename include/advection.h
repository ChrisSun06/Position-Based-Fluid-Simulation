#ifndef ADVECTION
#define ADVECTION
#include <Eigen/Core>
#include <Eigen/Sparse>
// Advection step of the fluid simulation.
//
// Inputs:
//   Position  The matrix containing position of fluid particles
//   Velocity  The matrix containing velocity of fluid particles
//   dt        Change in time
// Outputs:
//   Position  The matrix containing position of fluid particles
//
void advection(
  Eigen::MatrixXd & Position,
  const Eigen::MatrixXd & Velocity,
  const double dt);
#endif
