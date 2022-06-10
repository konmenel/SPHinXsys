#include "Dambreak_ch.h"

#define ENABLE_WATER 1

using std::flush;

int main()
{
	SPHSystem system(system_domain_bounds, resolution_ref);
	LogOutput fcout("Run.out");
	GlobalStaticVariables::physical_time_ = 0.0;

#if ENABLE_WATER
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
#if ENABLE_WATER
	ComplexBodyRelation water_block_complex(water_block, {&wall_boundary, &boulder});
#endif //ENABLE_WATER

	solid_dynamics::CorrectConfiguration boulder_corrected_configuration(boulder_inner);

	//-------------------------------------------------------------------
	//this section define all numerical methods will be used in this case
	//-------------------------------------------------------------------
	//-------- common particle dynamics ----------------------------------------
	Gravity gravity(Vec3d(0.0, 0.0, -gravity_g));
#if ENABLE_WATER
	TimeStepInitialization initialize_a_fluid_step(water_block, gravity);
	//-------- fluid dynamics --------------------------------------------------
	//evaluation of density by summation approach
	fluid_dynamics::DensitySummationFreeSurfaceComplex update_density_by_summation(water_block_complex);
	//time step size without considering sound wave speed
	fluid_dynamics::AdvectionTimeStepSize get_fluid_advection_time_step_size(water_block, U_f);
	//time step size with considering sound wave speed
	fluid_dynamics::AcousticTimeStepSize get_fluid_time_step_size(water_block);

	//pressure relaxation using verlet time stepping
	fluid_dynamics::PressureRelaxationWithWall pressure_relaxation(water_block_complex);
	fluid_dynamics::DensityRelaxationWithWall density_relaxation(water_block_complex);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAccelerationWithWall viscous_acceleration(water_block_complex);
	/** Impose transport velocity. */
	// fluid_dynamics::TransportVelocityCorrectionComplex transport_velocity_correction(water_block_complex);
	/** viscous acceleration and transport velocity correction can be combined because they are independent dynamics. */
	// CombinedInteractionDynamics viscous_acceleration_and_transport_correction(viscous_acceleration, transport_velocity_correction);
#endif //ENABLE_WATER

	fcout << "Setting up Chrono" << endl;
	ChSystemNSC ch_system;
	auto boulder_ch = addBoulderCh(ch_system);
	addWallsCh(ch_system);
	fcout << "Bodies added!" << endl;

	fcout << "Creating forces" << endl;
	// Set up the force and torque objects
	auto force_ch = chrono_types::make_shared<ChForce>();
	auto torque_ch = chrono_types::make_shared<ChForce>();
	force_ch->SetMode(ChForce::ForceType::FORCE);
	torque_ch->SetMode(ChForce::ForceType::TORQUE);
	fcout << "SetMode finished!" << endl;
	boulder_ch->AddForce(force_ch);
	boulder_ch->AddForce(torque_ch);
	fcout << "Forces added to the body!" << endl;
	force_ch->SetVrelpoint(ChVector<>(.0, .0, .0));

	ch_system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);
	ch_system.SetSolverType(ChSolver::Type::PSOR);
	ch_system.SetSolverMaxIterations(50);
	ch_system.Set_G_acc(ChVector<>(0.0, 0.0, -gravity_g));

	BoulderSystemForChrono boulder_for_chrono(boulder, "Boulder_chrono", boulder.body_shape_);
	ConstrainSolidBodyPartByChrono boulder_constain(boulder, boulder_for_chrono, boulder_ch);
	TotalForceOnSolidBodyPartForChrono force_on_boulder(boulder, boulder_for_chrono, boulder_ch, ch_system);

	fcout << "Chrono Setup Finished!" << endl;

	// Output system
	In_Output in_output(system);
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

	GlobalStaticVariables::physical_time_ = 0.0;
	system.restart_step_ = 0;
	size_t number_of_iterations = system.restart_step_;
	Real dt = 0.001;
	const Real end_time = 2.0;
	const Real out_dt = 0.01;
	const size_t report_steps = 100;
	const size_t restart_write_steps = 500;

	/** If the starting time is not zero, please setup the restart time step or read in restart states. */
	if (system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.restart_step_);
		restart_io.readFromFile(system.restart_step_);
		water_block.updateCellLinkedList();
		water_block_complex.updateConfiguration();
	}


	//statistics for computing time
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	fcout << "Main loop started..." << endl;
	while (GlobalStaticVariables::physical_time_ < end_time) {
		Real integration_time = 0.0;
		while (integration_time < out_dt) {
#if ENABLE_WATER
			// acceleration due to viscous force and gravity
			initialize_a_fluid_step.parallel_exec();
			Real Dt = get_fluid_advection_time_step_size.parallel_exec();
			Dt = SMIN(Dt, out_dt - integration_time);
			update_density_by_summation.parallel_exec();
			viscous_acceleration.parallel_exec();
			// transport_velocity_correction.parallel_exec();
			// viscous_acceleration_and_transport_correction.parallel_exec();
#else
			Real Dt = 0.005;
#endif // ENABLE_WATER

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) {
#if ENABLE_WATER
				pressure_relaxation.parallel_exec(dt);
				density_relaxation.parallel_exec(dt);
				dt = get_fluid_time_step_size.parallel_exec();
				dt = SMIN(dt, Dt - relaxation_time);
#endif // ENABLE_WATER
				SimTK::SpatialVec torque_force = force_on_boulder.parallel_exec();
				torque_ch->SetDir(vecToCh(torque_force[0]));
				force_ch->SetDir(vecToCh(torque_force[1]));
				torque_ch->SetMforce(torque_force[0].norm());
				force_ch->SetMforce(torque_force[1].norm());
				// boulder_ch->Accumulate_torque(vecToCh(torque_force[0]), false);
				// boulder_ch->Accumulate_force(vecToCh(torque_force[1]), ChVector<>(0.0, 0.0, 0.0), false);
				ch_system.DoStepDynamics(dt);
				boulder_constain.parallel_exec();
				
				integration_time += dt;
				relaxation_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}

			if (number_of_iterations % report_steps == 0) {
				fcout << "Step=" << number_of_iterations
				<< "\tTime=" << GlobalStaticVariables::physical_time_
				<< "\nforce_vec=" << force_ch->GetForce()
				<< "\tforce_modulus=" << force_ch->GetForceMod() << endl;
			}

			if (number_of_iterations % restart_write_steps == 0) {
				tick_count t2 = tick_count::now();
				restart_io.writeToFile(number_of_iterations);
				tick_count t3 = tick_count::now();
				interval += t3 - t2;
			}

			number_of_iterations++;

#if ENABLE_WATER
			water_block.updateCellLinkedList();
			water_block_complex.updateConfiguration();
#endif // ENABLE_WATER
		}

		tick_count t2 = tick_count::now();
		write_body_states.writeToFile(number_of_iterations);
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
