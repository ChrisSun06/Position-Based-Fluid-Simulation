#ifndef GRID_TO_PARTICLE
#define GRID_TO_PARTICLE
#include <Eigen/Core>
#include <grid_data.h>
// update Velocity based on updated grid velocity;
//
// Inputs:
//   grid      the grid;
//   Position  The matrix containing position of fluid particles
//   Velocity  The matrix containing velocity of fluid particles
// Outputs:
//   Velocity  The matrix containing velocity of fluid particles
//
void grid_to_particle(
    const grid_data & grid,
    const Eigen::MatrixXd & Position,
    Eigen::MatrixXd & Velocity);
#endif
