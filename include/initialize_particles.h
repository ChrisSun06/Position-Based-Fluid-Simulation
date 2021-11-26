#ifndef INITIALIZE_PARTICLES
#define INITIALIZE_PARTICLES
#include <Eigen/Core>
// initialize particle positions for simulation
//
// Inputs:
//   P: Matrix used to store point positions.
//   Starting_corner: starting corner for particles
//   step: distance between two neighbouting particles
//  nx: number of particles along X axis
//  ny: number of particles along Y axis
//  nz: number of particles along Z axis
//
void initialize_particles(
    Eigen::MatrixXd & P,
    const Eigen::RowVector3d Starting_corner,
    const double step,
    const int nx,
    const int ny,
    const int nz);
#endif
