#ifndef P_TO_GRID_INTERPOLATE
#define P_TO_GRID_INTERPOLATE
#include <Eigen/Core>
#include <Eigen/Sparse>
// Using trilinear interpolation to get stagger grid
//
// Inputs:
//   nx  number of grid steps along the x-direction
//   ny  number of grid steps along the y-direction
//   nz  number of grid steps along the z-direction
//   h  grid step size
//   corner  list of bottom-left-front corner position of grid
//   P  n by 3 list of query point locations
//   V  used to construct the stagger grid
//   direction: 0 x, 1 y, 2 z;
// Outputs:
//   Grid output grid;
//
void p_to_grid_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  const Eigen::VectorXd & V,
  const int direction,
  Eigen::VectorXd & Grid);
#endif
