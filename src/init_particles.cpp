#include <init_particles.h>
#include <iostream>
void init_particles(
	Particles& particles,
	Eigen::Vector3d& corner,
	Eigen::Vector3i& num_points,
	double step
) {
	
	int num_V = num_points(0) * num_points(1) * num_points(2);
	
	particles.position.resize(num_V, 3);
	particles.position.setZero();
	particles.PredictedPos.resize(num_V, 3);
	particles.PredictedPos.setZero();
	particles.deltaP.resize(num_V, 3);
	particles.deltaP.setZero();

	particles.velocity.resize(num_V, 3);
	particles.velocity.setZero();
	particles.acceleration.resize(num_V, 3);
	particles.acceleration.setZero();

	particles.pressure.resize(num_V);
	particles.pressure.setZero();
	particles.density.resize(num_V);
	particles.density.setZero();
	particles.force.resize(num_V, 3);
	particles.force.setZero();

	particles.C.resize(num_V);
	particles.C.setZero();
	particles.grad_C.resize(num_V);
	particles.grad_C.setZero();
	particles.accum_Grad_C.resize(num_V, 3);
	particles.accum_Grad_C.setZero();
	particles.lambda.resize(num_V);
	particles.lambda.setZero();

	for (int i = 0; i < num_points(0); i++) {
		for (int j = 0; j < num_points(1); j++) {
			for (int k = 0; k < num_points(2); k++) {
				particles.position.row(i * num_points(1) * num_points(2) + j * num_points(2) + k) <<
					corner(0) + step * i,
					corner(1) + step * j,
					corner(2) + step * k;
			}
		}
	}

	
}