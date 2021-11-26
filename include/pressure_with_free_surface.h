#ifndef PRESSURE_WITH_FREE_SURFACE
#define PRESSURE_WITH_FREE_SURFACE
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
void pressure_with_free_surface(
    grid_data & grid,
    const double rho,
    const double dt);
#endif
