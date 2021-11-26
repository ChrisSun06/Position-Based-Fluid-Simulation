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
	particles.PredictedPos.resize(num_V, 3);

	particles.velocity.resize(num_V, 3);
	particles.velocity.setZero();
	particles.acceleration.resize(num_V, 3);
	particles.acceleration.setZero();

	particles.pressure.resize(num_V);
	particles.density.resize(num_V);

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