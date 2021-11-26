#ifndef RAW_AIR_FLUID
#define RAW_AIR_FLUID
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <grid_data.h>
// marking cell as fluid cell or air cell
//
// Inputs:
//   P  n by 3 list of query point locations
//   grid grid data
// Outputs:
//   grid with marked air cells and fluid cells;
//
void raw_air_fluid(
  const Eigen::MatrixXd & P,
  grid_data & grid);
#endif
