#include "Dambreak_ch.h"

#define ENABLE_WATER


int main(int argc, char *argv[])
{
	SPHSystem system(system_domain_bounds, resolution_ref);
	LogOutput fcout("Run.out");

#ifdef BOOST_AVAILABLE
	system.handleCommandlineOptions(argc, argv);
#endif

	system.restart_step_ = 0;

	// Parameters for simulation
	GlobalStaticVariables::physical_time_ = 0.0;
	Real &sim_time = GlobalStaticVariables::physical_time_;
	size_t number_of_iterations = system.restart_step_;
	Real dt = 0.001;
	const Real end_time = 2.0;
	const Real out_dt = 0.01;
	const size_t report_steps = 100;
	const size_t min_restart_write_step = 500;	// Contact should occur after this step
	const size_t max_restart_write_step = 600;	// Contact should occure by this step

	// Creating the bodies.
#ifdef ENABLE_WATER
	//the water block
	WaterBlock water_block(system, "WaterBody");
	// create fluid particles
	FluidParticles fluid_particles(water_block, makeShared<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f));
#endif //ENABLE_WATER

	//the walls
	WallBoundary wall_boundary(system, "Wall");
	// create solid particles
	SolidParticles wall_particles(wall_boundary);

	//the boulder
	Boulder boulder(system, "Boulder");
	// create solid particles
	ElasticSolidParticles boulder_particles(boulder, makeShared<LinearElasticSolid>(rho0_s, poisson, Youngs_modulus));

	/** topology */
	BodyRelationInner boulder_inner(boulder);
#ifdef ENABLE_WATER
	ComplexBodyRelation water_block_complex(water_block, {&wall_boundary, &boulder});
	BodyRelationContact boulder_fluid_contact(boulder, { &water_block });
#endif //ENABLE_WATER

	solid_dynamics::CorrectConfiguration boulder_corrected_configuration(boulder_inner);

	//-------------------------------------------------------------------
	//this section define all numerical methods will be used in this case
	//-------------------------------------------------------------------
	//-------- common particle dynamics ----------------------------------------
	Gravity gravity(Vec3d(0.0, 0.0, -gravity_g));
#ifdef ENABLE_WATER
	TimeStepInitialization initialize_a_fluid_step(water_block, gravity);
	//-------- fluid dynamics --------------------------------------------------
	//evaluation of density by summation approach
	fluid_dynamics::DensitySummationFreeSurfaceComplex update_density_by_summation(water_block_complex);
	//time step size without considering sound wave speed
	fluid_dynamics::AdvectionTimeStepSize get_fluid_advection_time_step_size(water_block, U_f);
	//time step size with considering sound wave speed
	fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(water_block);

	//pressure relaxation using verlet time stepping
	fluid_dynamics::PressureRelaxationRiemannWithWall pressure_relaxation(water_block_complex);
	fluid_dynamics::DensityRelaxationRiemannWithWall density_relaxation(water_block_complex);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAccelerationWithWall viscous_acceleration(water_block_complex);
	/** Impose transport velocity. */
	// fluid_dynamics::TransportVelocityCorrectionComplex transport_velocity_correction(water_block_complex);
	/** viscous acceleration and transport velocity correction can be combined because they are independent dynamics. */
	// CombinedInteractionDynamics viscous_acceleration_and_transport_correction(viscous_acceleration, transport_velocity_correction);
	
	// Fluid-solid contact forces
	solid_dynamics::FluidPressureForceOnSolid fluid_pressure_force_on_boulder(boulder_fluid_contact);
	solid_dynamics::FluidViscousForceOnSolid fluid_viscous_force_on_boulder(boulder_fluid_contact);
#endif //ENABLE_WATER

	// Average velocity and acceleration on boulder
	// solid_dynamics::AverageVelocityAndAcceleration	average_velocity_and_acceleration(boulder);
	solid_dynamics::UpdateElasticNormalDirection 	boulder_update_normal(boulder);
	

	fcout << "Setting up Chrono" << endl;
	ChSystemNSC ch_system;
	auto boulder_ch = addBoulderCh(ch_system);
	addWallsCh(ch_system);
	fcout << "Bodies added!" << endl;

	ch_system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);
	ch_system.SetSolverType(ChSolver::Type::PSOR);
	ch_system.SetSolverMaxIterations(50);
	ch_system.Set_G_acc(ChVector<>(0.0, 0.0, -gravity_g));

	BoulderSystemForChrono boulder_for_chrono(boulder, "Boulder_chrono", boulder.body_shape_);
	ConstrainSolidBodyPartByChrono boulder_constain(boulder, boulder_for_chrono, boulder_ch);
	TotalForceOnSolidBodyPartForChrono force_on_boulder(boulder, boulder_for_chrono, boulder_ch, ch_system);

	fcout << "Chrono Setup Finished!" << endl;

	// Output system setup
	In_Output in_output(system);
	// Restart step should be set before the instanciation of In_Output because the
	// restart folder is deleted if restart set is 0!.
	RestartIO restart_io(in_output, system.real_bodies_);
	
	BodyStatesRecordingToLegacyVtk write_body_states(in_output, system.real_bodies_);
	write_body_states.writeToFile(0);

	/**
	 * @brief Setup geometrics and initial conditions
	 */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromBodyShape();
	boulder_particles.initializeNormalDirectionFromBodyShape();
	boulder_corrected_configuration.parallel_exec();

	/** If the starting time is not zero, please setup the restart time step or read in restart states. */
	if (system.restart_step_ != 0)
	{
		fcout << "Reading from restart files..." << endl;
		sim_time = restart_io.readRestartFiles(system.restart_step_);
			wall_boundary.updateCellLinkedList();
			boulder.updateCellLinkedList();
#ifdef ENABLE_WATER
			water_block.updateCellLinkedList();
			water_block_complex.updateConfiguration();
			boulder_fluid_contact.updateConfiguration();
#endif // ENABLE_WATER
		fcout << "Reading successful!\n"
			<< "Continuing from:\n"
			<< "Step=" << system.restart_step_
			<< "\tTime=" << sim_time << endl;
	}

	fcout << "Main loop started..." << endl;

	//statistics for computing time
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	while (sim_time < end_time) {
		Real integration_time = 0.0;
		while (integration_time < out_dt) {
#ifdef ENABLE_WATER
			// acceleration due to viscous force and gravity
			initialize_a_fluid_step.parallel_exec();
			Real Dt = get_fluid_advection_time_step_size.parallel_exec();
			Dt = SMIN(Dt, out_dt - integration_time);
			update_density_by_summation.parallel_exec();
			// viscous_acceleration.parallel_exec();
			// transport_velocity_correction.parallel_exec();
			// viscous_acceleration_and_transport_correction.parallel_exec();

			fluid_viscous_force_on_boulder.parallel_exec();
#else
			Real Dt = 0.005;
#endif // ENABLE_WATER
			boulder_update_normal.parallel_exec();

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) {
#ifdef ENABLE_WATER
				dt = get_fluid_time_step_size.parallel_exec();
				dt = SMIN(dt, Dt - relaxation_time);
				pressure_relaxation.parallel_exec(dt);
				fluid_pressure_force_on_boulder.parallel_exec();
				density_relaxation.parallel_exec(dt);
#endif // ENABLE_WATER

				// Solid dynamics
				// average_velocity_and_acceleration.initialize_displacement_.parallel_exec();
				
				SimTK::SpatialVec torque_force = force_on_boulder.parallel_exec();
				boulder_ch->Empty_forces_accumulators();
				boulder_ch->Accumulate_torque(vecToCh(torque_force[0]), false);
				boulder_ch->Accumulate_force(vecToCh(torque_force[1]), boulder_ch->GetPos(), false);
				ch_system.DoStepDynamics(dt);
				boulder_constain.parallel_exec();
				
				// average_velocity_and_acceleration.update_averages_.parallel_exec(dt);

				// Final step of timestepping, increment times
				integration_time += dt;
				relaxation_time += dt;
				sim_time += dt;
			}

			if (number_of_iterations % report_steps == 0) {
				fcout << "Step=" << number_of_iterations
				<< "\tTime=" << sim_time
				<< "\nForce_vec=" << boulder_ch->Get_accumulated_force() << endl;
			}

			wall_boundary.updateCellLinkedList();
			boulder.updateCellLinkedList();
#ifdef ENABLE_WATER
			water_block.updateCellLinkedList();
			water_block_complex.updateConfiguration();
			boulder_fluid_contact.updateConfiguration();
#endif // ENABLE_WATER
			
			number_of_iterations++;
		}

		tick_count t2 = tick_count::now();
		write_body_states.writeToFile(number_of_iterations);

		if (	number_of_iterations < max_restart_write_step
			&& 	number_of_iterations > min_restart_write_step)
		{
			restart_io.writeToFile(number_of_iterations);
		}

		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	fcout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
	fcout.close();

	return 0;
}
