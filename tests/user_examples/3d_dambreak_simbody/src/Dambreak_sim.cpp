/* ---------------------------------------------------------------------------*
*                       SPHinXsys: 3D dambreak example                        *
* ----------------------------------------------------------------------------*
* This is the one of the basic test cases for efficient and accurate time     *
* integration scheme investigation 							  				  *
* ---------------------------------------------------------------------------*/
#include "Dambreak_sim.h"

// TODO: Too slow, might be simbody contact.

//the main program
int main()
{

	//build up context -- a SPHSystem
	SPHSystem system(system_domain_bounds, resolution_ref);
	ofstream fcout("Run.out");

	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Tag for computation from restart files. 0: not from restart files. */
	system.restart_step_ = 0;

	//the water block
	WaterBlock water_block(system, "WaterBody");
	// create fluid particles
	FluidParticles fluid_particles(water_block, makeShared<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f));

	//the wall boundary
	WallBoundary wall_boundary(system, "Wall");
	//create solid particles
	SolidParticles wall_particles(wall_boundary);

	//the wall boundary
	Boulder boulder(system, "Boulder");
	//create solid particles
	ElasticSolidParticles boulder_particles(boulder, makeShared<LinearElasticSolid>(rho0_s, poisson, Youngs_modulus));
	
	/** topology */
	BodyRelationInner boulder_inner(boulder);
	ComplexBodyRelation water_block_complex(water_block, {&wall_boundary, &boulder});
	
	solid_dynamics::CorrectConfiguration boulder_corrected_configuration(boulder_inner);

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
	fluid_dynamics::ViscousAccelerationWithWall viscous_acceleration(water_block_complex);
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityCorrectionComplex transport_velocity_correction(water_block_complex);
	/** viscous acceleration and transport velocity correction can be combined because they are independent dynamics. */
	CombinedInteractionDynamics viscous_acceleration_and_transport_correction(viscous_acceleration, transport_velocity_correction);
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
		SimTK::Transform(Vec3d(VWx - BDx - BDL*0.5, BDy, BDz + BDH*0.5)),
		boulder_info, SimTK::Transform());

	SimTK::State state = MBsystem.realizeTopology();
	SimTK::RungeKutta2Integrator integ(MBsystem);
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

	size_t number_of_iterations = system.restart_step_;
	int screen_output_interval = 200;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 2.0;
	//time step size for output file
	Real D_Time = 0.01;
	Real dt = 0.001; //default acoustic time step sizes

	//statistics for computing time
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	fcout << "Start time stepping..." << endl;
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
			// transport_velocity_correction.parallel_exec();
			viscous_acceleration_and_transport_correction.parallel_exec();

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				pressure_relaxation.parallel_exec(dt);
				density_relaxation.parallel_exec(dt);
				dt = get_fluid_time_step_size.parallel_exec();
				dt = SMIN(dt, Dt - relaxation_time);

				SimTK::State& state_for_update = integ.updAdvancedState();
				force_on_bodies.clearAllBodyForces(state_for_update);
				force_on_bodies.setOneBodyForce(state_for_update, boulder_body, force_on_boulder.parallel_exec());
				integ.stepBy(dt);
				constraint_boulder.parallel_exec();

				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}

			if (number_of_iterations % screen_output_interval == 0)
			{
				fcout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
						  << GlobalStaticVariables::physical_time_
						  << "	dt = " << dt << endl;
			}
			
			number_of_iterations++;

			water_block.updateCellLinkedList();
			water_block_complex.updateConfiguration();
		}


		tick_count t2 = tick_count::now();
		write_water_block_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	fcout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
