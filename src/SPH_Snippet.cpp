#include "../includes/SPH_2D.h"
#include "../includes/file_writer.h"
#include <memory>
#include <string>
#include <sstream>
#include <iostream>
#include <omp.h>
#include <list>
#include <algorithm>

using namespace std;

shared_ptr<SPH_main> domain(new SPH_main()); // domain is an instance of the SPH_main class

//shared_ptr<SPH_main> updated_domain(new SPH_main()); // HOW to copy stuff inside domain into updated domain?

double SPH_main::h;

/**
choose dt to be the minimum of the three equation.
@param domain  the domain which contain the
*/
double find_optimal_time_step(shared_ptr<SPH_main> domain)
{
    list<double> mag_acc;
    list<double> all_rho;
    // loop through all the particles and find the particle with the lowest acceleration
    for (int i = 0; i < domain->particle_list.size(); i++)
    {
        // adding every particle's density value to the all_rho vector
        all_rho.push_back(domain->particle_list[i].rho);
        // finding the magnitude of every particle's acceleration and adding it to the mag_acc vector
        double acc = sqrt(pow(domain->particle_list[i].a[0], 2.0) + pow(domain->particle_list[i].a[1], 2.0));
        mag_acc.push_back(acc);
    }
    // finding the minimum density and acceleration
    list<double>::iterator max_particle_density = max_element(begin(all_rho), end(all_rho));
    list<double>::iterator max_particle_acc = max_element(begin(mag_acc), end(mag_acc));
    double max_rho = *max_particle_density;
    double max_acc = *max_particle_acc;
    double dt_F = sqrt(domain->dx/max_rho);
    double dt_A = domain->dx/(domain->c0*sqrt(pow(max_acc/1000, domain->gamma-1)));
    double C_CFL = 0.2; // can adjust this value (should be in range of 0.1-0.3)
    return C_CFL*min(dt_F, dt_A);
}

/**
implement the euler_time_step
@param domain the grid domain
@param dt the timestep
*/
void euler_time_step(shared_ptr<SPH_main> domain, double dt)
{
	// loop through every particle in the older domain:
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < domain->particle_list.size(); i++)
	{
		// only do euler with non boundary particles
		if (!domain->particle_list[i].boundary_particle)
		{
			// loop through x an y dimension
			for (int k = 0; k < 2; k++) // looping throgh vectors
			{
				// do simple euler step for position and velocity
				domain->particle_list[i].x[k] = domain->particle_list[i].x[k] + dt * (domain->particle_list[i].v[k]);
				domain->particle_list[i].v[k] = domain->particle_list[i].v[k] + dt * (domain->particle_list[i].a[k]);
			}
		}
		// euler step for density
		domain->particle_list[i].rho = domain->particle_list[i].rho + dt * (domain->particle_list[i].D);
	}
}
/**
implement the half time step
@param domian the gird domain
@param dt the timestep

*/
void half_full_time_step(shared_ptr<SPH_main> domain, double dt)
{
	// initiate new domain object to store half timestep results
	shared_ptr<SPH_main> half_domain(new SPH_main(*domain));

	// get half timestep results with simple euler
	euler_time_step(half_domain, dt/2);
	// update the grid 
	half_domain->allocate_to_grid();

	#pragma omp parallel
	{
		// do second step of improved euler with 
		#pragma omp parallel for schedule(static) no wait
		for (int i = 0; i < half_domain->particle_list.size(); i++)
			{
				half_domain->particle_list[i].rho = half_domain->particle_list[i].rho + dt * (half_domain->particle_list[i].D);
				if (isnan(half_domain->particle_list[i].rho))
					half_domain->particle_list[i].rho = 1000; 
				half_domain->particle_list[i].P = half_domain->particle_list[i].cal_pressure(half_domain->particle_list[i].rho, half_domain->gamma, half_domain->c0);

				half_domain->particle_list[i].mass = half_domain->particle_list[i].cal_mass(half_domain->dx, half_domain->particle_list[i].rho);
				half_domain->particle_list[i].calc_index();	
			}

		#pragma omp parallel for schedule(static)
		for (int i = 0; i < domain->particle_list.size(); i++)
		{
			half_domain->neighbour_iterate(&domain->particle_list[i]); 
			if (!domain->particle_list[i].boundary_particle)
			{
				for (int k = 0; k < 2; k++) // looping throgh vectors
				{
					domain->particle_list[i].x[k] = domain->particle_list[i].x[k] + dt * (half_domain->particle_list[i].v[k]);
					domain->particle_list[i].v[k] = domain->particle_list[i].v[k] + dt * (half_domain->particle_list[i].a[k]);
				}
			}
			domain->particle_list[i].rho = domain->particle_list[i].rho + dt * (domain->particle_list[i].D);
		}
	}
}

int main(void)
{

	system("exec rm results/*");

	double t = 0;

	int count = 0; // keeping track of the number of steps (iterations). needed because dt is changing between iterations

	int out_cnt = 0;

	while (t < domain->sim_time) // sim_length is the length of the simulation
	{
		if (domain->adaptive_timestepping)
			domain->dt = find_optimal_time_step(domain);
		domain->allocate_to_grid();

		if (count % 10 == 0)
		{
			// update rho function here
			for (int i = 0; i < domain->particle_list.size(); i++) {
				domain->particle_list[i].update_rho();
			}

			for (int i = 0; i < domain->particle_list.size(); i++)
			{
				domain->particle_list[i].rho = domain->particle_list[i].rho_new;
				domain->particle_list[i].P = domain->particle_list[i].cal_pressure(domain->particle_list[i].rho, domain->gamma, domain->c0);
			}
		}

		if (count % 10 == 0)
		{
			cout << "Wrting to file. Iteration: " << count << endl;
			stringstream filename;
			filename <<"./results/data_" << out_cnt << ".vtp";

			write_file(filename.str().c_str(), &domain->particle_list);
			out_cnt++;
			// this output needs to be the updated domain
		}

		for (int i = 0; i < domain->particle_list.size(); i++)
		{
			// domain->particle_list[i].rho = domain->particle_list[i].rho + dt * (domain->particle_list[i].D);
			if (isnan(domain->particle_list[i].rho))
				domain->particle_list[i].rho = 1000; 
			domain->particle_list[i].P = domain->particle_list[i].cal_pressure(domain->particle_list[i].rho, domain->gamma, domain->c0);

			domain->particle_list[i].mass = domain->particle_list[i].cal_mass(domain->dx, domain->particle_list[i].rho);
			domain->particle_list[i].calc_index();
		}

		for (int i = 0; i < domain->particle_list.size(); i++)
		{
			domain->neighbour_iterate(&domain->particle_list[i]);   
		}

	
		if (domain->time_step_meth == "Euler")
			euler_time_step(domain, domain->dt);
		else if (domain->time_step_meth == "Improved Euler")
			half_full_time_step(domain, domain->dt);
		else
			assert(1 == 0 && "not implemented yet. Valid time step methods are 'Euler' and 'Improved Euler' ");

		t += domain->dt; // increment the timer by dt

		count++;
	}

	return 0;
}
