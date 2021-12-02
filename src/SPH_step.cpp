#include <SPH_step.h>
#include <vector>
#include <weight_funcs.h>
#include <coefs.h>
#include <iostream>
#include <collision_response.h>

#define CONST_INV_REST_DENSITY 1.0
#define RELAXATION 0.01
#define PRESSURE_K 0.1
#define PRESSURE_N 6

Eigen::MatrixXd N(6, 3);
// N << 0., -1., 0.,
// 	0., 1., 0.,
// 	1., 0., 0.,
// 	0., 0., 1.,
// 	-1., 0., 0.,
// 	0., 0., -1.;
Eigen::MatrixXd P(6, 3);
// P <<
// 	0., 2., 0.,
// 	0., -1., 0.,
// 	-1., 0., 0.,
// 	0, 0., -1.,
// 	2., 0., 0.,
// 	0., 0., 2.;

double
	density,
	d2Cs, // Color surface Laplacian
	r; // positions distance
Eigen::Vector3d d,
	f_pressure, // pressure 
	f_visco, // viscosity 
	f_surface, // surface tension 
	f_g, // gravity 
	grad_press, // W_pressure gradient
	grad_poly, // W_poly gradient
	dCs; // Color surface normal 

double min_x;
double min_y;
double min_z;
double max_x;
double max_y;
double max_z;

int cell_x;
int cell_y;
int cell_z;

std::vector<std::vector<std::vector<std::vector<int> > > > grid(cell_x, std::vector<std::vector<std::vector<int> > >(cell_y, std::vector<std::vector<int> >(cell_z, std::vector<int>())));
int grid_x, grid_y, grid_z;


void addExtForce(double dt, const Eigen::Vector3d & extF, Particles& particles, Coef& coef){
	N << 0., -1., 0.,
	0., 1., 0.,
	1., 0., 0.,
	0., 0., 1.,
	-1., 0., 0.,
	0., 0., -1.;
	P <<
	0., 2., 0.,
	0., -1., 0.,
	-1., 0., 0.,
	0, 0., -1.,
	2., 0., 0.,
	0., 0., 2.;
	Eigen::Vector3d pi;
	Eigen::Vector3d vi;
	Eigen::Vector3d p_wall;
	Eigen::Vector3d N_wall;
	for (int i = 0; i < particles.position.rows(); i++) {
		for (int w = 0; w < P.rows(); w++) {
			pi = particles.PredictedPos.row(i);
			vi = particles.velocity.row(i);
			p_wall = P.row(w);
			N_wall = N.row(w);
			double d = (p_wall - pi).dot(N_wall);
			if (d > 0.) {
				particles.PredictedPos.row(i) += 2.0* d * N_wall;
				particles.velocity.row(i) -=  0.9 * vi.dot(N_wall) * N_wall;
			}
		}
	}
	for (int i = 0; i < particles.position.rows(); i++) {
		particles.velocity.row(i) += dt*(particles.force.row(i)+extF.transpose())/coef.MASS;
		particles.deltaP.row(i).setZero();
		particles.PredictedPos.row(i) = particles.position.row(i)+particles.velocity.row(i)*dt;
		particles.force.row(i).setZero();
	}
	// std::cout << "Predicted pos :\n" << particles.PredictedPos.row(0) << std::endl;
	// std::cout << "velo :\n" << particles.velocity.row(0) << std::endl;
	// std::cout << "extf :\n" << extF.transpose() << std::endl;
	for (int j = 0; j < particles.position.rows(); j++) {
		for (int k = 0; k < P.rows(); k++) {
			pi = particles.PredictedPos.row(j);
			vi = particles.velocity.row(j);
			p_wall = P.row(k);
			N_wall = N.row(k);
			double d = (p_wall - pi).dot(N_wall);
			if (d > 0.) {
				// printf("help\n");
				particles.PredictedPos.row(j) += 2.0* d * N_wall;
				particles.velocity.row(j) -=  0.9 * vi.dot(N_wall) * N_wall;
			}
		}
	}
}

void findNeighbors(std::vector< std::vector<int> > &neighbors){
	for (grid_x = 0; grid_x < cell_x; grid_x++) {
		for (grid_y = 0; grid_y < cell_y; grid_y++) {
			for (grid_z = 0; grid_z < cell_z; grid_z++) {
				assert(grid_x < grid.size());
				assert(grid_y < grid[grid_x].size());
				assert(grid_z < grid[grid_x ][grid_y ].size());
				for (int i : grid[grid_x][grid_y][grid_z]) {
					// neighbors.push_back(std::vector<int>());
					for (int count_x = -1; count_x < 2; count_x++) {
						if (grid_x + count_x < 0 || grid_x + count_x >= cell_x) continue;
						for (int count_y = -1; count_y < 2; count_y++) {
							if (grid_y + count_y < 0 || grid_y + count_y >= cell_y) continue;
							for (int count_z = -1; count_z < 2; count_z++) {
								if (grid_z + count_z < 0 || grid_z + count_z >= cell_z) continue;
								assert(grid_x + count_x >= 0 && grid_x + count_x < grid.size());
								assert(grid_y + count_y >= 0 && grid_y + count_y < grid[grid_x + count_x ].size());
								assert(grid_z + count_z >= 0 && grid_z + count_z < grid[grid_x + count_x][grid_y + count_y].size());
								for (int j : grid[grid_x + count_x][grid_y + count_y][grid_z + count_z]) {
									if(i != j){
										neighbors[i].push_back(j);
									}
								}
							}
						}
					}
				}
			}
		}
	}
}



void solveLambda(Particles& particles, Coef& coef, std::vector< std::vector<int> > neighbors){
	Eigen::Vector3d p;
	Eigen::Vector3d p_j;
	Eigen::Vector3d g;
	/////////////////Solve C_i///////////////////////////////////////////////
	for (int i = 0; i < particles.position.rows(); i++) {
		p = particles.PredictedPos.row(i);
		double ro=0.0;
		for(int j : neighbors[i]){
				if(i == j) continue;
				p_j = particles.PredictedPos.row(j);
				double kv=kernel(p-p_j,coef.H);
				particles.density(i) += kv * coef.MASS;
				particles.density(j) += kv * coef.MASS;
		}
		Eigen::Vector3d tmp;
		tmp.setZero();
		particles.density(i) += kernel(tmp,coef.H) * coef.MASS;
		// Particles[i]->density+=kernel(vec3(),core_radius)/Particles[i]->inv_mass;
	}

	for (int i = 0; i < particles.position.rows(); i++) {
		particles.C(i)=particles.density(i)*CONST_INV_REST_DENSITY-1.0;
		particles.density(i)=0.0;
	}
	/////////////////////////////////////////////////////////////////////////

	////////////////////////Solve C_i_gradient///////////////////////////////
	
	for (int i = 0; i < particles.position.rows(); i++) {
		p = particles.PredictedPos.row(i);
		double C_i_gradient, sum_gradients=0.0;
		for(int j : neighbors[i]){
				
			if(i == j) continue;
			p_j = particles.PredictedPos.row(j);
			spiky(p-p_j,coef.H, g);
			C_i_gradient=std::pow(g.norm()*CONST_INV_REST_DENSITY,2);
			//sum_gradients+=C_i_gradient;
			particles.grad_C(i)+=C_i_gradient;
			particles.grad_C(j)+=C_i_gradient;
		}

		for(int j : neighbors[i]){
			if(i == j) continue;
			p_j = particles.PredictedPos.row(j);
			Eigen::Vector3d agc;
			spiky(p-p_j,coef.H, agc);
			agc = agc * CONST_INV_REST_DENSITY;
			particles.accum_Grad_C.row(i)+=agc;
			particles.accum_Grad_C.row(j)-=agc;
		}
		C_i_gradient=particles.accum_Grad_C.row(i).norm()*CONST_INV_REST_DENSITY;
		particles.grad_C(i)+=C_i_gradient*C_i_gradient;

		double sum_Ci=particles.grad_C(i)+RELAXATION;
		particles.lambda(i)=-1.0*(particles.C(i)/sum_Ci);
		particles.grad_C(i)=0.0;
		particles.accum_Grad_C.row(i).setZero();
	}
}

void solvePosition(Particles& particles, Coef& coef, std::vector< std::vector<int> > neighbors){
	Eigen::Vector3d p;
	Eigen::Vector3d p_j;
	for (int i = 0; i < particles.position.rows(); i++) {
		p = particles.PredictedPos.row(i);
		double l=particles.lambda(i);
		Eigen::Vector3d delta;
		
		double k_term;
		Eigen::Vector3d tmp;
		tmp << 1.0, 1.0, 1.0;
		Eigen::Vector3d d_q= (0.1*coef.H)*tmp+p;

		double s_corr=0.0;

		for(int j : neighbors[i]){
				
			if(i == j) continue;
			p_j = particles.PredictedPos.row(j);
			double poly6pd_q=kernel(p-d_q,coef.H);

			if(poly6pd_q<0.1){ 
				k_term=0.0;
			} else {
				k_term= kernel(p-p_j,coef.H)/poly6pd_q;
			}
			s_corr = -PRESSURE_K * std::pow(k_term, PRESSURE_N);
			//s_corr=0.0f;
			//delta+=(l+(*pit)->lambda+s_corr)*spikyGradient(p-p_j,core_radius);
			Eigen::Vector3d agc;
			spiky(p-p_j,coef.H, agc);
			Eigen::Vector3d dp=(particles.lambda(i)+particles.lambda(j)+s_corr)* agc*CONST_INV_REST_DENSITY;
			particles.deltaP.row(i)+=dp;
			particles.deltaP.row(j)-=dp;
		}
	}
		

	for (int i = 0; i < particles.position.rows(); i++) {
		// if(Particles[i]->isOnSurface)
		// 	continue;
		particles.PredictedPos.row(i)+=particles.deltaP.row(i);
		particles.deltaP.row(i).setZero();
		/*
		In position based dynamics, collision handling is just moving particles back into a valid position, its
		velocity will be implied by the position change
		We add a wall to left, right and bottom. We leave top open for more fun. Potential problem: particles outside of
		the screen does not belong to any grid cell, therefore they don't have liquid properties.
		*/
		N << 0., -1., 0.,
		0., 1., 0.,
		1., 0., 0.,
		0., 0., 1.,
		-1., 0., 0.,
		0., 0., -1.;
		P <<
		0., 2., 0.,
		0., -1., 0.,
		-1., 0., 0.,
		0, 0., -1.,
		2., 0., 0.,
		0., 0., 2.;
		for (int w = 0; w < P.rows(); w++) {
			Eigen::Vector3d pi = particles.PredictedPos.row(i);
			Eigen::Vector3d vi = particles.velocity.row(i);
			Eigen::Vector3d p_wall = P.row(w);
			Eigen::Vector3d N_wall = N.row(w);
			double d = (p_wall - pi).dot(N_wall);

			if (d > 0.) {
				// push particle out of wall
				particles.PredictedPos.row(i) += 2.0 * d * N_wall;
				// reflect the velocity component # Don't need to do this!
				particles.velocity.row(i) -=  0.9 * vi.dot(N_wall) * N_wall;
			}
		}
	}
}

void final_update(Particles& particles, double dt){
	for (int i = 0; i < particles.position.rows(); i++) {
		// if(Particles[i]->isOnSurface)
		// 	continue;
		particles.velocity.row(i)=(particles.PredictedPos.row(i)-particles.position.row(i))/dt;
		particles.position.row(i)=particles.PredictedPos.row(i);
	}

	//apply_vorticity(dt);
	//apply_viscosity();
	// for(unsigned int i=0;i<ParticleCont;i++){
	// 	if(!dropToGround&&Particles[i]->isOnSurface)
	// 		continue;
	// 	Particles[i]->position=Particles[i]->PredictedPos;
	// }
}

void cleanup(Particles& particles){
	particles.PredictedPos.setZero();
	particles.deltaP.setZero();
	particles.C.setZero();
	particles.grad_C.setZero();
	particles.accum_Grad_C.setZero();
	particles.lambda.setZero();
	particles.density.setZero();
}

void PBF_update(Particles& particles, double dt, Coef& coef) {
	
	std::vector<int> v_row;
	std::vector< std::vector<int> > neighbors(particles.position.rows(), v_row);
	
	min_x = particles.position.col(0).minCoeff();
	min_y = particles.position.col(1).minCoeff();
	min_z = particles.position.col(2).minCoeff();

	printf("particles_min[0]=%f\n", min_x);
	printf("particles_min[1]=%f\n", min_y);
	printf("particles_min[2]=%f\n", min_z);

	f_g << 0, coef.G* 1.0, 0;

	max_x = particles.position.col(0).maxCoeff();
	max_y = particles.position.col(1).maxCoeff();
	max_z = particles.position.col(2).maxCoeff();

	printf("particles_max[0]=%f\n", max_x);
	printf("particles_max[1]=%f\n", max_y);
	printf("particles_max[2]=%f\n", max_z);

	cell_x = (max_x - min_x) / coef.H + 1;
	cell_y = (max_y - min_y) / coef.H + 1;
	cell_z = (max_z - min_z) / coef.H + 1;

	grid = std::vector<std::vector<std::vector<std::vector<int> > > >(cell_x, std::vector<std::vector<std::vector<int> > >(cell_y, std::vector<std::vector<int> >(cell_z, std::vector<int>())));
	cleanup(particles);

	// N <<
	// 0., -1., 0.,
	// 0., 1., 0.,
	// 1., 0., 0.,
	// 0., 0., 1.,
	// -1., 0., 0.,
	// 0., 0., -1.;
	// P <<
	// 0., 2., 0.,
	// 0., -1., 0.,
	// -1., 0., 0.,
	// 0, 0., -1.,
	// 2., 0., 0.,
	// 0., 0., 2.;

	// collision_response(particles, P, N);

	addExtForce(dt, f_g, particles, coef);
	for (size_t i = 0; i < particles.PredictedPos.rows(); i++) {
		Eigen::Vector3d pos = particles.PredictedPos.row(i);
		grid_x = (pos(0) - min_x) / coef.H;
		grid_y = (pos(1) - min_y) / coef.H;
		grid_z = (pos(2) - min_z) / coef.H;
		assert(grid_x < grid.size());
		assert(grid_y < grid[grid_x].size());
		assert(grid_z < grid[grid_x][grid_y].size());
		grid[grid_x][grid_y][grid_z].push_back(i);
	}

	findNeighbors(neighbors);

	for(unsigned int iter=0;iter<3;iter++){
		solveLambda(particles, coef, neighbors);
		solvePosition(particles, coef, neighbors);
	}

	final_update(particles, dt);
}

// void SPH_step(Particles& particles, double dt, Coef& coef) {


// 	// double
// 	// 	density,
// 	// 	d2Cs, // Color surface Laplacian
// 	// 	r; // positions distance
// 	// Eigen::Vector3d
// 	// 	d,
// 	// 	f_pressure, // pressure 
// 	// 	f_visco, // viscosity 
// 	// 	f_surface, // surface tension 
// 	// 	f_g, // gravity 
// 	// 	grad_press, // W_pressure gradient
// 	// 	grad_poly, // W_poly gradient
// 	// 	dCs; // Color surface normal 

// 	min_x = particles.position.col(0).minCoeff();
// 	min_y = particles.position.col(1).minCoeff();
// 	min_z = particles.position.col(2).minCoeff();

// 	max_x = particles.position.col(0).maxCoeff();
// 	max_y = particles.position.col(1).maxCoeff();
// 	max_z = particles.position.col(2).maxCoeff();

// 	cell_x = (max_x - min_x) / coef.H + 1;
// 	cell_y = (max_y - min_y) / coef.H + 1;
// 	cell_z = (max_z - min_z) / coef.H + 1;

// 	for (size_t i = 0; i < particles.position.rows(); i++) {
// 		Eigen::Vector3d pos = particles.position.row(i);
// 		grid_x = (pos(0) - min_x) / coef.H;
// 		grid_y = (pos(1) - min_y) / coef.H;
// 		grid_z = (pos(2) - min_z) / coef.H;
// 		grid[grid_x][grid_y][grid_z].push_back(i);
// 	}

// 	for (grid_x = 0; grid_x < cell_x; grid_x++) {
// 		for (grid_y = 0; grid_y < cell_y; grid_y++) {
// 			for (grid_z = 0; grid_z < cell_z; grid_z++) {
// 				for (int i : grid[grid_x][grid_y][grid_z]) {
// 					density = 0;
// 					for (int count_x = -1; count_x < 2; count_x++) {

// 						if (grid_x + count_x < 0 || grid_x + count_x >= cell_x) continue;

// 						for (int count_y = -1; count_y < 2; count_y++) {

// 							if (grid_y + count_y < 0 || grid_y + count_y >= cell_y) continue;

// 							for (int count_z = -1; count_z < 2; count_z++) {

// 								if (grid_z + count_z < 0 || grid_z + count_z >= cell_z) continue;

// 								for (int j : grid[grid_x + count_x][grid_y + count_y][grid_z + count_z]) {
// 									d = particles.position.row(i) - particles.position.row(j);
// 									r = sqrt(d.dot(d));
// 									if (r > coef.H) continue;
// 									density += Wpoly(r, coef.H);
// 								}
// 							}
// 						}
// 					}

// 					particles.density(i) = density * coef.MASS;
// 					particles.pressure(i) = coef.STIFFNESS * (particles.density(i) - coef.RHO_IDEAL);
// 				}
// 			}
// 		}
// 	}


// 	for (grid_x = 0; grid_x < cell_x; grid_x++) {
// 		for (grid_y = 0; grid_y < cell_y; grid_y++) {
// 			for (grid_z = 0; grid_z < cell_z; grid_z++) {

// 				for (int i : grid[grid_x][grid_y][grid_z]) {

// 					// forces and surface tension

// 					f_pressure.setZero();
// 					f_surface.setZero();
// 					f_g << 0, coef.G* particles.density(i), 0;
// 					f_visco.setZero();
// 					dCs.setZero();
// 					d2Cs = 0;

// 					for (int count_x = -1; count_x < 2; count_x++) {

// 						if (grid_x + count_x < 0 || grid_x + count_x >= cell_x) continue;

// 						for (int count_y = -1; count_y < 2; count_y++) {

// 							if (grid_y + count_y < 0 || grid_y + count_y >= cell_y) continue;

// 							for (int count_z = -1; count_z < 2; count_z++) {

// 								if (grid_z + count_z < 0 || grid_z + count_z >= cell_z) continue;

// 								for (int j : grid[grid_x + count_x][grid_y + count_y][grid_z + count_z]) {
// 									d = particles.position.row(i) - particles.position.row(j);
// 									r = sqrt(d.dot(d));
// 									if (r > coef.H) continue;
// 									if (r > 0.) {
// 										dWpress(d, r, coef.H, grad_press);
// 										f_pressure -= coef.MASS * (particles.pressure(i) + particles.pressure(j)) / (2. * particles.density(j)) * grad_press;
// 										Eigen::Vector3d diff_v = particles.velocity.row(j) - particles.velocity.row(i);
// 										f_visco += d2Wvisco(r, coef.H) * diff_v / particles.density(j);
// 									}
// 									dWpoly(d, r, coef.H, grad_poly);
// 									dCs += particles.density(j) * grad_poly;
// 									d2Cs += d2Wpoly(r, coef.H) / particles.density(j);
// 								}
// 							}
// 						}
// 					}

// 					f_visco *= coef.VISCOSITY * coef.MASS;
// 					dCs *= coef.MASS;
// 					d2Cs *= coef.MASS;



// 					if (dCs.norm() > coef.SUFRACE_THRESH) {
// 						f_surface = -coef.SUFRACE_TENSION * d2Cs * dCs.normalized();
// 					}
// 					particles.acceleration.row(i) = (f_pressure + f_visco + f_surface + f_g) / particles.density(i);
// 					particles.velocity.row(i) += dt * particles.acceleration.row(i);
// 					particles.position.row(i) += dt * particles.velocity.row(i);
// 				}
// 			}
// 		}
// 	}
// }
