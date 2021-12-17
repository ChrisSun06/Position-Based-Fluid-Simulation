#ifndef KERNELS
#define KERNELS

#include <Eigen/Dense>

double kernel(Eigen::Ref<const Eigen::Vector3d> r, double H);

void spiky(
	Eigen::Ref<const Eigen::Vector3d> r, double H, Eigen::Vector3d& g);

double Wpoly(double r, double H);

#endif