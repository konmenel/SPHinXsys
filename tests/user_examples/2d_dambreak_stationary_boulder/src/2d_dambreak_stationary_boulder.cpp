/**
 * @file 	Dambreak_boulder.cpp
 * @brief 	2D dambreak example with moving obsacle.
 * @details This is a case to help me understand how to code works.
 * @author 	Constantinos Menelaou
*/
/**
 * @brief The setup of the case
 */
#include "case.h"
/**
 * @brief 	Main program starts here.
 */


int main()
{
	/** Build up context -- a SPHSystem. */
	SPHSystem system(system_domain_bounds, particle_spacing_ref);

	/** Define external force.*/
	Gravity gravity(Vecd(0.0, -gravity_g));
	/** The water block, body, material and particles container. */
	WaterBlock *water_block 		= new WaterBlock(system, "WaterBody");
	WaterMaterial *water_material 	= new WaterMaterial();
	FluidParticles fluid_particles(water_block, water_material);
	/** The wall boundary, body and particles container. */
	WallBoundary *wall_boundary 	= new WallBoundary(system, "Wall");
	SolidParticles wall_particles(wall_boundary);
	/** Boulder system. Body, material and particle container. */
	Boulder *boulder 					= new Boulder(system, "Boulder");
	SolidParticles boulder_particles(boulder);
	/** topology */
	ComplexBodyRelation* water_block_complex = new ComplexBodyRelation(water_block, { wall_boundary, boulder });
	ContactBodyRelation* boulder_fluid_contact 		= new ContactBodyRelation(boulder, { water_block });
	/** 
	 * Methods used for time stepping
	 */
	/** Time step initialization, add gravity. */
	InitializeATimeStep 	initialize_gravity_to_fluid(water_block, &gravity);
	/** Evaluation of density by summation approach. */
	fluid_dynamics::DensitySummationFreeSurfaceComplex update_density_by_summation(water_block_complex);
	/** time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize	get_fluid_advection_time_step_size(water_block, U_f);
	/** time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize		get_fluid_time_step_size(water_block);
	/** pressure relaxation using verlet time stepping. */
	fluid_dynamics::PressureRelaxationRiemannWithWall	pressure_relaxation(water_block_complex);
	fluid_dynamics::DensityRelaxationRiemannWithWall	density_relaxation(water_block_complex);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAccelerationWithWall viscous_acceleration(water_block_complex);
	
	/** Output. */
	In_Output in_output(system);
	WriteBodyStatesToVtu 		write_real_body_states(in_output, system.real_bodies_);
	WriteTotalForceOnSolid      write_total_force_on_boulder(in_output, boulder);
	/**
	 * @brief Prepare quantities will be used once only and initial condition.
	 */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromGeometry();
	boulder_particles.initializeNormalDirectionFromGeometry();

	write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
	write_total_force_on_boulder.WriteToFile(GlobalStaticVariables::physical_time_);
	/** Simulation start here. */
	/** starting time zero. */
	GlobalStaticVariables::physical_time_ = 0.0;
	system.restart_step_ = 0;
	int number_of_iterations = 0;
	int screen_output_interval = 500;
	Real End_Time = 2.0;					/**< End time. */
	Real D_Time = End_Time / 1100.0;		/**< Time stamps for output of body states. */
	Real Dt = 0.0;							/**< Default advection time step sizes. */
	Real dt = 0.0; 							/**< Default acoustic time step sizes. */
	Real total_time = 0.0;
	/** statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/** Main Loop. */
	write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
	while (GlobalStaticVariables::physical_time_ < End_Time){
		// std::cout << "Iteration " << number_of_iterations << "\n";
		Real integral_time = 0.0;
		while (integral_time < D_Time) {
			/** Acceleration due to viscous force and gravity. */
			initialize_gravity_to_fluid.parallel_exec();
			
			Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_density_by_summation.parallel_exec();
			viscous_acceleration.parallel_exec();
			
			/** Dynamics including pressure relaxation. */
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) {
				pressure_relaxation.parallel_exec(dt);
				density_relaxation.parallel_exec(dt);
				
				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integral_time += dt;
				total_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			
			if (number_of_iterations % screen_output_interval == 0) {
				std::cout << fixed << setprecision(9) << "N=" << number_of_iterations 
					 << "	Total Time = " << total_time 
					 << "	Physical Time = " << GlobalStaticVariables::physical_time_ 
					 << "	Dt = " << Dt << "	dt = " << dt << "\n";
			}
			number_of_iterations++;
			// damping_wave.parallel_exec(Dt);
			water_block->updateCellLinkedList();
			wall_boundary->updateCellLinkedList();
			boulder->updateCellLinkedList();
			water_block_complex->updateConfiguration();
			boulder_fluid_contact->updateConfiguration();
			write_total_force_on_boulder.WriteToFile(GlobalStaticVariables::physical_time_);
		}

		tick_count t2 = tick_count::now();
		write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_ / (D_Time * 1e4));
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds.\n";
	
	return 0;
}
