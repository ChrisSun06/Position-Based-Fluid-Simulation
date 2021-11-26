#include "advection.h"

void advection(
    Eigen::MatrixXd & Position,
    const Eigen::MatrixXd & Velocity,
    const double dt)
{
    Position = Position + dt * Velocity;
}
