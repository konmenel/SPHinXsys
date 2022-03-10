/**
 * @file 	Dambreak_boulder.cpp
 * @brief 	2D dambreak example with moving obsacle.
 * @details This is a case to help me understand how to code works.
 * @author 	Constantinos Menelaou
 */
/**
 * @brief The setup of the case
 */
#include "fine.h"
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
	WaterBlock water_block(system, "WaterBody");
	FluidParticles fluid_particles(water_block, makeShared<WaterMaterial>());
	/** The wall boundary, body and particles container. */
	WallBoundary wall_boundary(system, "Wall");
	SolidParticles wall_particles(wall_boundary);

	ObserverBody observer(system, "BoulderObserver");
	ObserverParticles observer_particles(observer, makeShared<ObserverParticleGenerator>());
	/** topology */
	ComplexBodyRelation water_block_complex(water_block, {&wall_boundary});
	BodyRelationContact observer_contact_with_water(observer, {&water_block});
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

	/** Output. */
	ofstream fcout("./stdout.out");
	fcout << "Setting up output...\n";
	In_Output in_output(system);
	fluid_particles.addAVariableToWrite<indexScalar, Real>("Pressure");
	BodyStatesRecordingToVtp write_real_body_states(in_output, system.real_bodies_);
	ObservedQuantityRecording<indexScalar, Real> pressure_probe("Pressure", in_output, observer_contact_with_water);
	/**
	 * @brief Prepare quantities will be used once only and initial condition.
	 */
	/** offset particle position */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();

	write_real_body_states.writeToFile(0);
	pressure_probe.writeToFile(0);
	/** Simulation start here. */
	fcout << "Setting up main loop...";
	/** starting time zero. */
	system.restart_step_ = 0;
	GlobalStaticVariables::physical_time_ = 0.0;
	int number_of_iterations = 0;
	int screen_output_interval = 1000;
	Real End_Time = 10.0;			 	/**< End time. */
	Real D_Time = 0.002; 				/**< Time stamps for output of body states. */
	Real Dt = 0.0;					 	/**< Default advection time step sizes. */
	Real dt = 0.0;					 	/**< Default acoustic time step sizes. */
	Real total_time = 0.0;
	/** statistics for computing time. */
	fcout << " DONE\n"
		  << "Main loop started.\n";
	fcout.flush();
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/** Main Loop. */
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integral_time = 0.0;
		while (integral_time < D_Time)
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
				density_relaxation.parallel_exec(dt);

				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integral_time += dt;
				total_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				fcout << fixed << setprecision(9) << "N=" << number_of_iterations
					 << "	Total Time = " << total_time
					 << "	Physical Time = " << GlobalStaticVariables::physical_time_
					 << "	Dt = " << Dt << "	dt = " << dt << "\n";
				fcout.flush();
			}
			number_of_iterations++;
			water_block.updateCellLinkedList();
			wall_boundary.updateCellLinkedList();
			water_block_complex.updateConfiguration();

			tick_count t2 = tick_count::now();
			pressure_probe.writeToFile(GlobalStaticVariables::physical_time_);
			tick_count t3 = tick_count::now();
			interval += t3 - t2;
		}

		tick_count t2 = tick_count::now();
		write_real_body_states.writeToFile(number_of_iterations);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t wall_time, total_run_time;
	wall_time = t4 - t1 - interval;
	total_run_time = t4 - t1;
	cout << "Total wall time for computation: " << wall_time.seconds() << " seconds.\n";
	cout << "Total time: " << total_run_time.seconds() << "seconds.\n";

	fcout.flush();
	fcout.close();
	return 0;
}
