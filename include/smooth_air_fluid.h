#ifndef SMOOTH_AIR_FLUID
#define SMOOTH_AIR_FLUID
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <grid_data.h>
// smooth Raw_air_fluid using 3d gaussian kernel
//
// Inputs:
//   grid grid data
// Outputs:
//   grid with smoothed air cells and fluid cells;
//
void smooth_air_fluid(
  grid_data & grid);
#endif
