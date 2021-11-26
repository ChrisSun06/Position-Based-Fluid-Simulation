#ifndef SOLID_BOUNDARY
#define SOLID_BOUNDARY
#include <Eigen/Core>
#include <grid_data.h>
// take care of velocity change on the boundary
//
// Inputs:
//   grid      the output grid;
// Outputs:
//   grid      the output grid;
//
void solid_boundary(
    grid_data & grid);
#endif
