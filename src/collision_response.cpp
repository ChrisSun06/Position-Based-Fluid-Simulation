#include <collision_response.h>
#include <iostream>

void collision_response(
	Particles &particles,
	Eigen::Ref<const Eigen::MatrixXd> P,
	Eigen::Ref<const Eigen::MatrixXd> N
) {
	for (int i = 0; i < particles.position.rows(); i++) {
		for (int w = 0; w < P.rows(); w++) {
			Eigen::Vector3d pi = particles.position.row(i);
			Eigen::Vector3d vi = particles.velocity.row(i);
			Eigen::Vector3d p_wall = P.row(w);
			Eigen::Vector3d N_wall = N.row(w);
			double d = (p_wall - pi).dot(N_wall);

			if (d > 0.) {
				// push particle out of wall
				particles.position.row(i) += d * N_wall;
				// reflect the velocity component
				particles.velocity.row(i) -=  vi.dot(N_wall) * N_wall;
			}
		}
	}

}