#ifndef PRESSURE
#define PRESSURE
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <grid_data.h>
// Calculate pressure;
//
// Inputs:
//   grid      regular and stagger grid;
//   rho       fluid density
//   dt        change in time
// Outputs:
//   grid      the output grid, with calculated pressure;
//
void pressure(
    grid_data & grid,
    const double rho,
    const double dt);
#endif
