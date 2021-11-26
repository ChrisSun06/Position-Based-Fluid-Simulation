#ifndef PRESSURE_PROJECTION
#define PRESSURE_PROJECTION
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <grid_data.h>
// Use pressure to update velocity, after this step X, Y, and Z of the grid would store the change in velocity;
//
// Inputs:
//   grid      regular and stagger grid;
//   rho       fluid density
//   dt        change in time
// Outputs:
//   grid      the output grid, with updated velocity;
//
void pressure_projection(
    grid_data & grid,
    const double rho,
    const double dt);
#endif
