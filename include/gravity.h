#ifndef GRAVITY
#define GRAVITY
#include <Eigen/Core>
#include <Eigen/Sparse>
// external force step of the fluid simulation.
//
// Inputs:
//   Velocity  The matrix containing velocity of fluid particles
//   g         g
//   dt        Change in time
// Outputs:
//   Velocity  The matrix containing velocity of fluid particles
//
void gravity(
  Eigen::MatrixXd & Velocity,
  const Eigen::Vector3d & g,
  const double dt);
#endif
