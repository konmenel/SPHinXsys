/* ---------------------------------------------------------------------------*
*                       SPHinXsys: 3D dambreak example                        *
* ----------------------------------------------------------------------------*
* This is the one of the basic test cases for efficient and accurate time     *
* integration scheme investigation 							  				  *
* ---------------------------------------------------------------------------*/
#include "Dambreak_sim.h"

#define ENABLE_WATER 1

// TODO: Too slow, might be simbody contact.

//the main program
int main(int argc, char *argv[])
{

	//build up context -- a SPHSystem
	SPHSystem system(system_domain_bounds, resolution_ref);
	LogOutput fcout("Run.out");

#ifdef BOOST_AVAILABLE
	system.handleCommandlineOptions(argc, argv);
#endif

	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for computation from restart files. 0: not from restart files. */
	system.restart_step_ = 0;

	Real &sim_time = GlobalStaticVariables::physical_time_;
	size_t number_of_iterations = system.restart_step_;
	Real dt = 0.001;
	const Real end_time = 2.0;
	const Real out_dt = 0.01;
	const size_t report_steps = 100;
	const size_t min_restart_write_step = 500;	// Contact should occur after this step
	const size_t max_restart_write_step = 600;	// Contact should occure by this step

#if ENABLE_WATER
	//the water block
	WaterBlock water_block(system, "WaterBody");
	// create fluid particles
	FluidParticles fluid_particles(water_block, makeShared<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f));
#endif //ENABLE_WATER

	//the wall boundary
	WallBoundary wall_boundary(system, "Wall");
	//create solid particles
	SolidParticles wall_particles(wall_boundary);
	fcout << "Tank created!\n";

	//the wall boundary
	fcout << "Creating boulder..." << endl;
	Boulder boulder(system, "Boulder");
	//create solid particles
	ElasticSolidParticles boulder_particles(boulder, makeShared<LinearElasticSolid>(rho0_s, poisson, Youngs_modulus));
	fcout << "Boulder created!\n";
	
	/** topology */
	BodyRelationInner boulder_inner(boulder);
#if ENABLE_WATER
	ComplexBodyRelation water_block_complex(water_block, {&wall_boundary, &boulder});
	BodyRelationContact boulder_fluid_contact(boulder, { &water_block });
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
	solid_dynamics::AverageVelocityAndAcceleration	average_velocity_and_acceleration(boulder);
	solid_dynamics::UpdateElasticNormalDirection 	boulder_update_normal(boulder);
	

	fcout << "Simbody starting..." << endl;
	// ----------------------------------------------------------------------------
	//Simbody
	//-----------------------------------------------------------------------------
	/** Multi-body system. */
	/** set up the multi body system. */
	SimTK::MultibodySystem MBsystem;
	/** the bodies or matter of the system. */
	SimTK::SimbodyMatterSubsystem 	matter(MBsystem);
	/** the forces of the system. */
	SimTK::GeneralForceSubsystem 	forces(MBsystem);
	/** the contact syster. */
	SimTK::ContactTrackerSubsystem  tracker(MBsystem);
	SimTK::CompliantContactSubsystem contactForces(MBsystem, tracker);
	contactForces.setTransitionVelocity(1e-3);

	//-----------------------------------------------------------------------------
	//Simbody forces
	//-----------------------------------------------------------------------------
	SimTK::Force::UniformGravity sim_gravity(forces, matter, Vec3d(0.0, 0.0, -gravity_g));
	/** discrete forces acting on the bodies. */
	SimTK::Force::DiscreteForces force_on_bodies(forces, matter);
	
	/** contanct clique for wall */
	SimTK::ContactCliqueId clique = SimTK::ContactSurface::createNewContactClique();
	
	addSimbodyWallContacts(matter, clique);
	addCliffContactForSimbody(matter, clique);
	
	BoulderSystemForSimbody boulder_multibody(boulder, "BoulderMultibody", boulder.body_shape_);
	SimTK::Body::Rigid    boulder_info(*boulder_multibody.body_part_mass_properties_);

	addBoulderContactForSimbody(boulder_info);

	SimTK::MobilizedBody::Free boulder_body(matter.Ground(), 
		Vec3d(VWx - BDx - 0.5*BDL, BDy, BDz + 0.5*BDH),
		boulder_info, SimTK::Transform());

	SimTK::State state = MBsystem.realizeTopology();
	SimTK::SemiExplicitEuler2Integrator integ(MBsystem);
	integ.setAccuracy(1e-3);
	integ.setAllowInterpolation(false);
	integ.initialize(state);

	solid_dynamics::TotalForceOnSolidBodyPartForSimBody
		force_on_boulder(boulder, boulder_multibody, MBsystem, boulder_body, force_on_bodies, integ);
	solid_dynamics::ConstrainSolidBodyPartBySimBody
		constraint_boulder(boulder, boulder_multibody, MBsystem, boulder_body, force_on_bodies, integ);

	//-----------------------------------------------------------------------------
	//outputs
	//-----------------------------------------------------------------------------
	In_Output in_output(system);
	BodyStatesRecordingToLegacyVtk write_water_block_states(in_output, system.real_bodies_);
	/** Output the body states for restart simulation. */
	RestartIO restart_io(in_output, system.real_bodies_);

	//-------------------------------------------------------------------
	//from here the time stepping begins
	//-------------------------------------------------------------------

	/**
	 * @brief Setup geometrics and initial conditions
	 */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromBodyShape();
	boulder_particles.initializeNormalDirectionFromBodyShape();
	boulder_corrected_configuration.parallel_exec();
	/**
	* @brief The time stepping starts here.
	*/
	/** If the starting time is not zero, please setup the restart time step or read in restart states. */
	// if (system.restart_step_ != 0)
	// {
	// 	GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(system.restart_step_);
	// 	water_block.updateCellLinkedList();
	// 	water_block_complex.updateConfiguration();
	// }

	/** Output the start states of bodies. */
	write_water_block_states.writeToFile(0);

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

	//statistics for computing time
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	fcout << "Start time stepping..." << endl;
	//computation loop starts
	while (sim_time < end_time)
	{
		Real integration_time = 0.0;
		//integrate time (loop) until the next output time
		while (integration_time < out_dt)
		{

#if ENABLE_WATER
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
			while (relaxation_time < Dt)
			{
#if ENABLE_WATER
				dt = get_fluid_time_step_size.parallel_exec();
				dt = SMIN(dt, Dt - relaxation_time);
				
				pressure_relaxation.parallel_exec(dt);
				fluid_pressure_force_on_boulder.parallel_exec();
				density_relaxation.parallel_exec(dt);
#else
				dt = 0.0002;
#endif // ENABLE_WATER

				// Solid dynamics
				average_velocity_and_acceleration.initialize_displacement_.parallel_exec();
				
				SimTK::State& state_for_update = integ.updAdvancedState();
				force_on_bodies.clearAllBodyForces(state_for_update);
				force_on_bodies.setOneBodyForce(state_for_update, boulder_body, force_on_boulder.parallel_exec());
				integ.stepBy(dt);
				constraint_boulder.parallel_exec();
				
				average_velocity_and_acceleration.update_averages_.parallel_exec(dt);

				// Final step of timestepping, increment times
				relaxation_time += dt;
				integration_time += dt;
				sim_time += dt;
			}

			if (number_of_iterations % report_steps == 0)
			{
				fcout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << sim_time
						  << "	dt = " << dt << endl;
			}
			
			wall_boundary.updateCellLinkedList();
			boulder.updateCellLinkedList();
#if ENABLE_WATER
			water_block.updateCellLinkedList();
			water_block_complex.updateConfiguration();
			boulder_fluid_contact.updateConfiguration();
#endif // ENABLE_WATER

			number_of_iterations++;
		}


		tick_count t2 = tick_count::now();
		write_water_block_states.writeToFile(number_of_iterations);

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
