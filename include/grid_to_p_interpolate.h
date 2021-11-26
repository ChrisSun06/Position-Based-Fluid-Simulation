#ifndef GRID_TO_P_INTERPOLATE
#define GRID_TO_P_INTERPOLATE
#include <Eigen/Core>
#include <Eigen/Sparse>
// Use trilinear interpolation to update V, final step of pressure projection.
//
// Inputs:
//   nx  number of grid steps along the x-direction
//   ny  number of grid steps along the y-direction
//   nz  number of grid steps along the z-direction
//   h  grid step size
//   corner  list of bottom-left-front corner position of grid
//   P  n by 3 list of query point locations
//   V  to be updated values
//   direction: 0 x, 1 y, 2 z;
// Outputs:
//   V  updated values;
//
void grid_to_p_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::VectorXd & V,
  const int direction,
  const Eigen::VectorXd & Grid);
#endif
