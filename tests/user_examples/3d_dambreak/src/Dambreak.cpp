/* ---------------------------------------------------------------------------*
*                       SPHinXsys: 3D dambreak example                        *
* ----------------------------------------------------------------------------*
* This is the one of the basic test cases for efficient and accurate time     *
* integration scheme investigation 							  				  *
* ---------------------------------------------------------------------------*/
#include "Dambreak.h"


//the main program
int main()
{

	//build up context -- a SPHSystem
	SPHSystem system(system_domain_bounds, resolution_ref);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for computation from restart files. 0: not from restart files. */
	system.restart_step_ = 0;

	//the water block
	WaterBlock water_block(system, "WaterBody");
	//create fluid particles
	FluidParticles fluid_particles(water_block, makeShared<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f));

	//the wall boundary
	WallBoundary wall_boundary(system, "Wall");
	//create solid particles
	SolidParticles wall_particles(wall_boundary);

	//the wall boundary
	Building building(system, "Building");
	//create solid particles
	SolidParticles building_particles(building);

	ObserverBody fluid_observer(system, "Fluidobserver");
	//create observer particles
	ObserverParticles observer_particles(fluid_observer, makeShared<ObserverParticleGenerator>());

	/** topology */
	ComplexBodyRelation water_block_complex(water_block, {&wall_boundary, &building});
	BodyRelationContact fluid_observer_contact(fluid_observer, {&water_block});

	//-------------------------------------------------------------------
	//this section define all numerical methods will be used in this case
	//-------------------------------------------------------------------
	//-------- common particle dynamics ----------------------------------------
	Gravity gravity(Vec3d(0.0, 0.0, -gravity_g));
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
	// fluid_dynamics::ViscousAccelerationWithWall viscous_acceleration(water_block_complex);
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityCorrectionComplex transport_velocity_correction(water_block_complex);
	/** viscous acceleration and transport velocity correction can be combined because they are independent dynamics. */
	// CombinedInteractionDynamics viscous_acceleration_and_transport_correction(viscous_acceleration, transport_velocity_correction);

	//-----------------------------------------------------------------------------
	//outputs
	//-----------------------------------------------------------------------------
	In_Output in_output(system);
	BodyStatesRecordingToVtp write_water_block_states(in_output, system.real_bodies_);
	/** Output the body states for restart simulation. */
	RestartIO restart_io(in_output, system.real_bodies_);
	/** Output the mechanical energy of fluid body. */
	BodyReducedQuantityRecording<TotalMechanicalEnergy>
		write_water_mechanical_energy(in_output, water_block, gravity);
	ObservedQuantityRecording<indexScalar, Real>
		write_recorded_water_pressure("Pressure", in_output, fluid_observer_contact);
	//-------------------------------------------------------------------
	//from here the time stepping begins
	//-------------------------------------------------------------------

	/**
	 * @brief Setup geometrics and initial conditions
	 */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromBodyShape();
	/**
	* @brief The time stepping starts here.
	*/
	/** If the starting time is not zero, please setup the restart time step or read in restart states. */
	if (system.restart_step_ != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.restart_step_);
		water_block.updateCellLinkedList();
		water_block_complex.updateConfiguration();
	}

	/** Output the start states of bodies. */
	write_water_block_states.writeToFile(0);
	/** Output the mechanical energy of fluid. */
	write_water_mechanical_energy.writeToFile(0);

	size_t number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 1.6;
	//time step size for output file
	Real D_Time = 0.01;
	Real dt = 0.0; //default acoustic time step sizes

	//statistics for computing time
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	//computation loop starts
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		//integrate time (loop) until the next output time
		while (integration_time < D_Time)
		{

			//acceleration due to viscous force and gravity
			initialize_a_fluid_step.parallel_exec();
			Real Dt = get_fluid_advection_time_step_size.parallel_exec();
			Dt = SMIN(Dt, D_Time - integration_time);
			update_density_by_summation.parallel_exec();
			// viscous_acceleration.parallel_exec();
			// viscous_acceleration_and_transport_correction.parallel_exec();
			transport_velocity_correction.parallel_exec();

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				pressure_relaxation.parallel_exec(dt);
				density_relaxation.parallel_exec(dt);
				dt = get_fluid_time_step_size.parallel_exec();
				dt = SMIN(dt, Dt - relaxation_time);
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	Dt = " << Dt << "	dt = " << dt << "\n";

				if (number_of_iterations % restart_output_interval == 0)
					restart_io.writeToFile(number_of_iterations);
			}
			number_of_iterations++;

			water_block.updateCellLinkedList();
			water_block_complex.updateConfiguration();
			fluid_observer_contact.updateConfiguration();
			write_recorded_water_pressure.writeToFile(number_of_iterations);
		}

		write_water_mechanical_energy.writeToFile(number_of_iterations);

		tick_count t2 = tick_count::now();
		write_water_block_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
