#include <iostream>
#include "../includes/SPH_2D.h"
#include "../inlcudes/file_writer.h"

using namespace std;



void euler_time_step(shared_ptr<SPH_main> domain, double dt)
{
	// loop through every particle in the older domain:
#pragma omp parallel for num_threads(4)
	for (int i = 0; i < domain->particle_list.size(); i++)
	{
		if (!domain->particle_list[i].boundary_particle)
		{
			for (int k = 0; k < 2; k++) // looping throgh vectors
			{
				domain->particle_list[i].x[k] = domain->particle_list[i].x[k] + dt * (domain->particle_list[i].v[k]);
				domain->particle_list[i].v[k] = domain->particle_list[i].v[k] + dt * (domain->particle_list[i].a[k]);
			}
			// find value of x vector one step forward:
		   // x at new time-step = x at previous time-step + ( dt * velocity at previous time-step )

		   // find value of v vector one step forward
		   // v at new time-step = v at previous time-step + ( dt * acceleration at previous time-step )

		}
		// find value of density one step forward
		// rho at new time-step = rho at previous time-step + ( dt * D at previous time-step )
		domain->particle_list[i].rho = domain->particle_list[i].rho + dt * (domain->particle_list[i].D);
		domain->particle_list[i].P = domain->particle_list[i].cal_pressure(domain->particle_list[i].rho, domain->gamma, domain->c0);
		domain->particle_list[i].mass = domain->particle_list[i].cal_mass(domain->dx, domain->particle_list[i].rho);
		domain->particle_list[i].calc_index();
	}
}



void half_full_time_step(shared_ptr<SPH_main> domain, shared_ptr<SPH_main> updated_domain, double dt)
{
	
	euler_time_step(domain, half_domain, dt/2);
	particle_deflecter(half_domain);
	half_domain->allocate_to_grid();
	for (int i = 0; i < domain->particle_list.size(); i++)
	{
		domain->particle_list[i].mass = domain->particle_list[i].cal_mass(domain->dx, updated_domain->particle_list[i].rho);

		domain->particle_list[i].P = domain->particle_list[i].cal_pressure(domain->rho_0, domain->B, domain->gamma, domain->c0);

		domain->particle_list[i].calc_index();
		domain->neighbour_iterate(&domain->particle_list[i]);

		domain->particle_list[i].D = domain->particle_list[i].navier_rho(domain->particle_list[i].neighbours, domain->dx, SPH_main::dW_dr);

		domain->particle_list[i].navier_acc(domain->particle_list[i].neighbours, domain->dx, SPH_main::dW_dr, domain->gamma, domain->c0, domain->mu);
	}

	for (int i = 0; i < domain->particle_list.size(); i++)
    {
		if (!domain->particle_list[i].boundary_particle)
		{
			for (int k = 0; k < 2; k++) // looping throgh vectors
			{
				updated_domain->particle_list[i].x[k] = domain->particle_list[i].x[k] + dt/2 * (half_domain->particle_list[i].v[k]);
				updated_domain->particle_list[i].v[k] = domain->particle_list[i].v[k] + dt/2 * (half_domain->particle_list[i].a[k]);
			}
		}
		updated_domain->particle_list[i].rho = domain->particle_list[i].rho + dt/2 * (half_domain->particle_list[i].D);
		
		// updated_domain->particle_list[i].P = updated_domain->particle_list[i].cal_pressure(updated_domain->particle_list[i].rho, updated_domain->gamma, updated_domain->c0);
		
		// updated_domain->particle_list[i].mass = updated_domain->particle_list[i].cal_mass(updated_domain->dx, updated_domain->particle_list[i].rho);

		// updated_domain->particle_list[i].calc_index();
		// updated_domain->neighbour_iterate(&domain->particle_list[i]);
    }
}


void paticle_delfector(shared_ptr<SPH_main> updated_domain, double dt)
{
    for (int i=0; i < updated_domain->particle_list.size(); i++)
    {
        if (updated_domain->particle_list[i].x[0] > updated_domain->min_x[0] + updated_domain->dx)
        {//if go near the left wall
            updated_domain->particle_list[i].v[0] += dt*20;
        }
            
        if (updated_domain->particle_list[i].x[0] > updated_domain->max_x[0] - updated_domain->dx)
        {//if go near the right wall
            updated_domain->particle_list[i].v[0] -= dt*20;
        }
        
        if (updated_domain->particle_list[i].x[1] < updated_domain->min_x[1] + updated_domain->dx)
        {//if go near the bottom
            updated_domain->particle_list[i].v[1] += dt*20;
        }
            
        if (updated_domain->particle_list[i].x[1] > updated_domain->max_x[1] - updated_domain->dx)
        {//if go near the top
            updated_domain->particle_list[i].v[1] -= dt*20;
        }
            
            //deal with the particles that goes out of the boundry
        if (updated_domain->particle_list[i].x[0] < 0)
        {
            updated_domain->particle_list[i].x[0] = -updated_domain->particle_list[i].x[0]/2;
            for(int j=0;j<2;j++)
            {
            updated_domain->particle_list[i].v[j] = 0;
            }
        }
        
        if (updated_domain->particle_list[i].x[1] < 0)
        {
            updated_domain->particle_list[i].x[1] = 0;
            for(int j=0;j<2;j++)
            {
            updated_domain->particle_list[i].v[j] = 0;
            }
        }
        
        if (updated_domain->particle_list[i].x[0] > updated_domain->max_x[0])
        {
            updated_domain->particle_list[i].x[0] = updated_domain->max_x[0] - updated_domain->particle_list[i].x[0] + updated_domain->max_x[0];
            for(int j=0;j<2;j++)
            {
            updated_domain->particle_list[i].v[j] = 0;
            }
        }
        
        if (updated_domain->particle_list[i].x[1] > updated_domain->max_x[1])
        {
            updated_domain->particle_list[i].x[1] = updated_domain->max_x[1] - updated_domain->particle_list[i].x[1] + updated_domain->max_x[1];
            for(int j=0;j<2;j++)
            {
            updated_domain->particle_list[i].v[j] = 0;
            }
        }
    }
}









// void update(shared_ptr<SPH_main> updated_domain)
// {
// 	for (int i = 0; i < updated_domain->particle_list.size(); i++)
// 	{
// 		if (!updated_domain->particle_list[i].boundary_particle)
// 		{
// 			//deal with the particles that goes near the boundry
// 			if (updated_domain->particle_list[i].x[0] < 0.1) { //if go near the left wall
// 				updated_domain->particle_list[i].v[0] += dt * 20;
// 			}
// 			if (updated_domain->particle_list[i].x[0] > 20 - 0.1) {//if go near the right wall
// 				updated_domain->particle_list[i].v[0] -= dt * 20;
// 			}
// 			if (updated_domain->particle_list[i].x[1] < 0.1) {//if go near the bottom
// 				updated_domain->particle_list[i].v[1] += dt * 20;
// 			}
// 			if (updated_domain->particle_list[i].x[1] > 10 - 0.1) {//if go near the top
// 				updated_domain->particle_list[i].v[1] -= dt * 20;
// 			}

// 			//deal with the particles that goes out of the boundry
// 			if (updated_domain->particle_list[i].x[0] < 0) {
// 				updated_domain->particle_list[i].x[0] = -updated_domain->particle_list[i].x[0] / 2;
// 				for (int j = 0; j < 2; j++) {
// 					updated_domain->particle_list[i].v[j] = 0;
// 				}
// 			}
// 			if (updated_domain->particle_list[i].x[1] < 0) {
// 				updated_domain->particle_list[i].x[1] = 0;
// 				for (int j = 0; j < 2; j++) {
// 					updated_domain->particle_list[i].v[j] = 0;
// 				}
// 			}
// 			if (updated_domain->particle_list[i].x[0] > 20)
// 			{
// 				updated_domain->particle_list[i].x[0] = 20 - updated_domain->particle_list[i].x[0] + 20;
// 				for (int j = 0; j < 2; j++)
// 					updated_domain->particle_list[i].v[j] = 0;
// 			}
// 			if (updated_domain->particle_list[i].x[1] > 10)
// 			{
// 				updated_domain->particle_list[i].x[1] = 10 - updated_domain->particle_list[i].x[1] + 10;
// 				for (int j = 0; j < 2; j++)
// 					updated_domain->particle_list[i].v[j] = 0;
// 			}
// 		}
// 	}
// }




