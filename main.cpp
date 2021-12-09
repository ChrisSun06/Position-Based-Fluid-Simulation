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

int main(int argc, char* argv[])
{

	Eigen::Vector3i num_points(15, 15, 10);
	//Eigen::Vector3i num_points(10, 10, 10);
	double step_size = 0.05;
	Eigen::Vector3d corner(-0.5, -0.5, -0.5);
	const Eigen::RowVector3d particle_color(0.2, 0.8, 0.8);

	igl::opengl::glfw::Viewer viewer;
	const int xid = viewer.selected_data_index;
	viewer.append_mesh();
	Particles particles;
	Coef coef;
	coef.MASS = 1.0;
	coef.H = 0.09;
	coef.G = -9.8;


	init_particles(particles, corner, num_points, step_size);
	const auto& reset = [&]()
	{
		init_particles(particles, corner, num_points, step_size);
		viewer.data_list[xid].set_points(particles.position, (1. - (1. - particle_color.array()) * .9));
	};

	const auto& move = [&]()
	{
		// collision_response(particles, P, N);
		// SPH_step(particles, 0.02, coef);
		PBF_update(particles, 0.02, coef);

		viewer.data_list[xid].set_points(particles.position, (1. - (1. - particle_color.array()) * .9));

	};


	viewer.callback_key_pressed =
		[&](igl::opengl::glfw::Viewer&, unsigned char key, int)->bool
	{
		if(key == 'r' || key == 'R'){
			reset();
		}
		return false;
	};
	viewer.core().is_animating = true;
	viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer & )->bool
	{
		move();
		viewer.data_list[xid].set_points(particles.position, (1. - (1. - particle_color.array()) * .9));
		return false;
	};
	viewer.data_list[xid].set_colors(particle_color);
	viewer.data_list[xid].set_points(particles.position, (1. - (1. - particle_color.array()) * .9));
	viewer.data_list[xid].point_size = 6.0;
	viewer.launch();




	return EXIT_SUCCESS;
}


