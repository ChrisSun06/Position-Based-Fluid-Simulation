#include <particle_sys_data.h>
#include <init_particles.h>
#include <SPH_step.h>
#include <coefs.h>

#include <igl/list_to_matrix.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <unistd.h>

#define X_POINTS 15
#define Y_POINTS 15
#define Z_POINTS 15

Eigen::Vector3d start_point;
Particles particles;

void initialize_particles(Particles& particles){
	particles.position.resize(X_POINTS*Y_POINTS*Z_POINTS, 3);
	particles.position.setZero();
	particles.PredictedPos.resize(X_POINTS*Y_POINTS*Z_POINTS, 3);
	particles.PredictedPos.setZero();
	particles.deltaP.resize(X_POINTS*Y_POINTS*Z_POINTS, 3);
	particles.deltaP.setZero();

	particles.velocity.resize(X_POINTS*Y_POINTS*Z_POINTS, 3);
	particles.velocity.setZero();
	particles.acceleration.resize(X_POINTS*Y_POINTS*Z_POINTS, 3);
	particles.acceleration.setZero();

	particles.pressure.resize(X_POINTS*Y_POINTS*Z_POINTS);
	particles.pressure.setZero();
	particles.density.resize(X_POINTS*Y_POINTS*Z_POINTS);
	particles.density.setZero();
	particles.force.resize(X_POINTS*Y_POINTS*Z_POINTS, 3);
	particles.force.setZero();

	particles.C.resize(X_POINTS*Y_POINTS*Z_POINTS);
	particles.C.setZero();
	particles.grad_C.resize(X_POINTS*Y_POINTS*Z_POINTS);
	particles.grad_C.setZero();
	particles.accum_Grad_C.resize(X_POINTS*Y_POINTS*Z_POINTS, 3);
	particles.accum_Grad_C.setZero();
	particles.lambda.resize(X_POINTS*Y_POINTS*Z_POINTS);
	particles.lambda.setZero();
}

void reset_particles(Particles& particles){
	int count = 0;
	double step = 0.05;
	initialize_particles(particles);
	for (int x=0; x<X_POINTS; x++){
		for (int y=0; y<Y_POINTS; y++){
			for (int z=0; z<Z_POINTS; z++){
				particles.position.row(count) << x*step+start_point(0), y*step+start_point(1), z*step+start_point(2);
				count += 1;
			}
		}
	}
}

void update_particles(Particles& particles){
	PBF_update(particles, 0.03);
}

int main(int argc, char* argv[])
{

	const Eigen::RowVector3d particle_color(0.0, 0.3, 0.6);

	start_point << -0.5, -0.5, -0.5;

	igl::opengl::glfw::Viewer viewer;
	viewer.append_mesh();

	initialize_particles(particles);
	reset_particles(particles);
	viewer.data_list[0].set_points(particles.position, particle_color);

	viewer.callback_key_pressed =
		[&](igl::opengl::glfw::Viewer&, unsigned char key, int)->bool
	{
		if(key == 'r' || key == 'R'){
			reset_particles(particles);
			viewer.data_list[0].set_points(particles.position, particle_color);
		}
		return false;
	};
	viewer.core().is_animating = true;
	viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer & )->bool
	{
		update_particles(particles);
		viewer.data_list[0].set_points(particles.position, particle_color);
		return false;
	};
	viewer.data_list[0].set_colors(particle_color);
	viewer.data_list[0].set_points(particles.position, particle_color);
	viewer.data_list[0].point_size = 8.0;
	viewer.launch();




	return 0;
}


