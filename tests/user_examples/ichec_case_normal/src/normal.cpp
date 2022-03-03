/**
 * @file 	Dambreak_boulder.cpp
 * @brief 	2D dambreak example with moving obsacle.
 * @details This is a case to help me understand how to code works.
 * @author 	Constantinos Menelaou
 */
/**
 * @brief The setup of the case
 */
#include "coarse.h"
/**
 * @brief 	Main program starts here.
 */

// TODO: fix bouncing/found out why bounces higher
// TODO: fix contact dumping
// TODO: add dumping area at the end
// TODO: add observers, wave probes etc.

int main()
{
	/** Build up context -- a SPHSystem. */
	SPHSystem system(system_domain_bounds, particle_spacing_ref);

	/** Define external force.*/
	Gravity gravity(Vecd(0.0, -gravity_g));
	/** The water block, body, material and particles container. */
	WaterBlock water_block(system, "WaterBody");
	FluidParticles fluid_particles(water_block, makeShared<WaterMaterial>());
	/** The wall boundary, body and particles container. */
	WallBoundary wall_boundary(system, "Wall");
	// SolidParticles wall_particles(wall_boundary);
	/** topology */
	ComplexBodyRelation water_block_complex(water_block, {&wall_boundary});
	/**
	 * Methods used for time stepping
	 */
	/** Time step initialization, add gravity. */
	TimeStepInitialization initialize_gravity_to_fluid(water_block, gravity);
	/** Evaluation of density by summation approach. */
	fluid_dynamics::DensitySummationFreeSurfaceComplex update_density_by_summation(water_block_complex);
	/** time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize get_fluid_advection_time_step_size(water_block, U_f);
	/** time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(water_block);
	/** pressure relaxation using verlet time stepping. */
	fluid_dynamics::PressureRelaxationRiemannWithWall pressure_relaxation(water_block_complex);
	fluid_dynamics::DensityRelaxationRiemannWithWall density_relaxation(water_block_complex);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAccelerationWithWall viscous_acceleration(water_block_complex);

	/**
	 * @brief Prepare quantities will be used once only and initial condition.
	 */
	/** offset particle position */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();

	/** Simulation start here. */
	/** starting time zero. */
	system.restart_step_ = 0;
	GlobalStaticVariables::physical_time_ = 0.0;
	int number_of_iterations = 0;
	int screen_output_interval = 1000;
	Real End_Time = 5.0;			 /**< End time. */
	Real Dt = 0.0;					 /**< Default advection time step sizes. */
	Real dt = 0.0;					 /**< Default acoustic time step sizes. */
	Real total_time = 0.0;
	/** statistics for computing time. */
	tick_count t1 = tick_count::now();
	/** Main Loop. */
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		/** Acceleration due to viscous force and gravity. */
		initialize_gravity_to_fluid.parallel_exec();

		Dt = get_fluid_advection_time_step_size.parallel_exec();
		update_density_by_summation.parallel_exec();
		viscous_acceleration.parallel_exec();

		/** Dynamics including pressure relaxation. */
		Real relaxation_time = 0.0;
		while (relaxation_time < Dt)
		{
			pressure_relaxation.parallel_exec(dt);
			// boulder_damping.parallel_exec(dt);
			density_relaxation.parallel_exec(dt);

			dt = get_fluid_time_step_size.parallel_exec();
			relaxation_time += dt;
			total_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
		}

		if (number_of_iterations % screen_output_interval == 0)
		{
			cout << fixed << setprecision(9) << "N=" << number_of_iterations
					<< "	Total Time = " << total_time
					<< "	Physical Time = " << GlobalStaticVariables::physical_time_
					<< "	Dt = " << Dt << "	dt = " << dt << "\n";
		}
		number_of_iterations++;
		water_block.updateCellLinkedList();
		wall_boundary.updateCellLinkedList();
		water_block_complex.updateConfiguration();
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t wall_time, total_run_time;
	wall_time = t4 - t1;
	total_run_time = t4 - t1;
	cout << "Total wall time for computation: " << wall_time.seconds() << " seconds.\n";
	cout << "Total time: " << total_run_time.seconds() << "seconds.\n";
	return 0;
}
