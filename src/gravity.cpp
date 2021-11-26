#include "gravity.h"

void gravity(
  Eigen::MatrixXd & Velocity,
  const Eigen::Vector3d & g,
  const double dt)
{
    Eigen::VectorXd ones = Eigen::VectorXd::Ones(Velocity.rows());
    Eigen::MatrixXd G = ones * g.transpose();
    Velocity = Velocity + dt * G;
}
