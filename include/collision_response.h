#ifndef COLLISION_H
#define COLLISION_H

#include <Eigen/Dense>
#include <particle_sys_data.h>


void collision_response(
	Particles& particles,
	Eigen::Ref<const Eigen::MatrixXd> P,
	Eigen::Ref<const Eigen::MatrixXd> N
);

#endif // !COLLISION_H
