#ifndef PARTICLE_TO_GRID
#define PARTICLE_TO_GRID
#include <Eigen/Core>
#include <grid_data.h>
// Construct stagger grid;
//
// Inputs:
//   Position  The matrix containing position of fluid particles
//   Velocity  The matrix containing velocity of fluid particles
//   nx
//   ny
//   nz
//   h         the size of the grid cell
// Outputs:
//   grid      the output grid;
//
void particle_to_grid(
    grid_data & grid,
    const Eigen::MatrixXd & Position,
    const Eigen::MatrixXd & Velocity,
    const Eigen::RowVector3d & corner,
    const int nx,
    const int ny,
    const int nz,
    const double h);
#endif
