#include <Particles.h>
#include <PBF_solver.h>

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
Particles data;

void initialize_data(Particles& data){
	data.q.resize(X_POINTS*Y_POINTS*Z_POINTS, 3);
	data.q.setZero();
	data.q_new.resize(X_POINTS*Y_POINTS*Z_POINTS, 3);
	data.q_new.setZero();
	data.deltaP.resize(X_POINTS*Y_POINTS*Z_POINTS, 3);
	data.deltaP.setZero();

	data.v.resize(X_POINTS*Y_POINTS*Z_POINTS, 3);
	data.v.setZero();

	data.density.resize(X_POINTS*Y_POINTS*Z_POINTS);
	data.density.setZero();
	data.f.resize(X_POINTS*Y_POINTS*Z_POINTS, 3);
	data.f.setZero();

	data.C.resize(X_POINTS*Y_POINTS*Z_POINTS);
	data.C.setZero();
	data.grad_C.resize(X_POINTS*Y_POINTS*Z_POINTS);
	data.grad_C.setZero();
	data.lambda.resize(X_POINTS*Y_POINTS*Z_POINTS);
	data.lambda.setZero();
}

void reset_data(Particles& data){
	int count = 0;
	double step = 0.05;
	initialize_data(data);
	for (int x=0; x<X_POINTS; x++){
		for (int y=0; y<Y_POINTS; y++){
			for (int z=0; z<Z_POINTS; z++){
				data.q.row(count) << x*step+start_point(0), y*step+start_point(1), z*step+start_point(2);
				count += 1;
			}
		}
	}
}

void update_data(Particles& data){
	PBF_update(data, 0.03);
}

int main(int argc, char** argv){

	const Eigen::RowVector3d color(0.4, 0.7, 1);

	start_point << -0.5, -0.5, -0.5;

	igl::opengl::glfw::Viewer viewer;
	viewer.append_mesh();

	initialize_data(data);
	reset_data(data);
	viewer.data_list[0].set_points(data.q, color);

	viewer.callback_key_pressed =
		[&](igl::opengl::glfw::Viewer&, unsigned char key, int)->bool
	{
		if(key == 'r' || key == 'R'){
			reset_data(data);
			viewer.data_list[0].set_points(data.q, color);
		}
		return false;
	};
	viewer.core().is_animating = true;
	viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer & )->bool
	{
		update_data(data);
		viewer.data_list[0].set_points(data.q, color);
		return false;
	};
	viewer.data_list[0].set_points(data.q, color);
	viewer.data_list[0].point_size = 7.0;
	viewer.launch();




	return 0;
}


