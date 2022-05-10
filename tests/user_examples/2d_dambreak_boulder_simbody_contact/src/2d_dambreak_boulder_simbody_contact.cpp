/**
 * @file 	Dambreak_boulder.cpp
 * @brief 	2D dambreak example with moving obsacle.
 * @details This is a case to help me understand how to code works.
 * @author 	Constantinos Menelaou
*/
/**
 * @brief The setup of the case
 */
#include "2d_dambreak_boulder.h"
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
	FluidParticles fluid_particles(water_block, makeShared<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f));
	/** The wall boundary, body and particles container. */
	WallBoundary wall_boundary(system, "Wall");
	SolidParticles wall_particles(wall_boundary);
	/** Boulder system. Body, material and particle container. */
	Boulder boulder(system, "Boulder");
	ElasticSolidParticles boulder_particles(boulder, makeShared<LinearElasticSolid>(rho0_s, poisson, Youngs_modulus));
	/** topology */
	BodyRelationInner boulder_inner(boulder);
	ComplexBodyRelation water_block_complex(water_block, { &wall_boundary, &boulder });
	BodyRelationContact boulder_fluid_contact(boulder, { &water_block });
	/**
	 * Methods only used only once
	 */
	/** corrected strong configuration. */
	solid_dynamics::CorrectConfiguration boulder_corrected_configuration(boulder_inner);
	/** 
	 * Methods used for time stepping
	 */
	/** Time step initialization, add gravity. */
	TimeStepInitialization initialize_gravity_to_fluid(water_block, gravity);
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
	/** Inflow boundary condition. */
	// fluid_dynamics::DampingBoundaryCondition	damping_wave(water_block, new DampingBuffer(water_block, "DampingBuffer"));
	/** Fluid force on boulder. */
	solid_dynamics::FluidPressureForceOnSolid fluid_pressure_force_on_boulder(boulder_fluid_contact);
	/** Fluid viscous force on boulder. */
	solid_dynamics::FluidViscousForceOnSolid fluid_viscous_force_on_boulder(boulder_fluid_contact);
	/** average velocity for boulder. */
	solid_dynamics::AverageVelocityAndAcceleration	average_velocity_and_acceleration(boulder);
	solid_dynamics::UpdateElasticNormalDirection 	boulder_update_normal(boulder);
	
	/** Multi-body system. */
	cout << "Simbody setup... ";
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
	cout << "SUCCESS!\n";

	cout << "Adding forces... ";
	/** 
	 * @details Add gravity to mb body. 
	 * @param[in,out] forces, The subsystem to which this force should be added.
	 * @param[in]     matter, The subsystem containing the bodies that will be affected.
	 * @param[in]    gravity, The default gravity vector v, interpreted as v=g*d where g=|\a gravity| is 
     *				a positive scalar and d is the "down" direction unit vector d=\a gravity/g.
	 * @param[in]  zeroHeight This is an optional specification of the default value for the height
     *				up the gravity vector that is considered to be "zero" for purposes of
     *				calculating the gravitational potential energy. The default is 
     *				\a zeroHeight == 0, i.e., a body's potential energy is defined to be zero
     *				when the height of its mass center is the same as the height of the Ground
     *				origin. The zero height will have the value specified here unless
     *				explicitly changed within a particular State use the setZeroHeight() 
     *			method. 
	 * @par Force Each body B that has not been explicitly excluded will experience a force 
	 *		fb = mb*g*d, applied to its center of mass, where mb is the mass of body B. 
	 * @par Potential Energy
	 *		Gravitational potential energy for a body B is mb*g*hb where hb is the height of 
	 *		body B's mass center over an arbitrary "zero" height hz (default is hz=0), 
	 *		measured along the "up" direction -d. If pb is the Ground frame vector giving 
	 *		the position of body B's mass center, its height over or under hz is 
	 *		hb=pb*(-d) - hz. Note that this is a signed quantity so the potential energy is 
	 *		also signed. 0.475
	 */
	SimTK::Force::UniformGravity sim_gravity(forces, matter, Vec3d(0.0, -gravity_g, 0.0));
	/** discrete forces acting on the bodies. */
	SimTK::Force::DiscreteForces force_on_bodies(forces, matter);
	cout << "SUCCESS!\n";
	
	/** Wall contacts */
	cout << "Adding wall contacts... ";
	/** contanct clique for wall */
	SimTK::ContactCliqueId clique = SimTK::ContactSurface::createNewContactClique();
	
	addSimbodyWallContacts(matter, clique);

	// Cliff body
	addCliffContactForSimbody(matter, clique);
	cout << "SUCCESS!\n";

	/** mass properties of boulder. */
	cout << "Adding boulder... ";
	MultiPolygonShape flap_multibody_shape(createBoulderSimbodyConstrainShape());
	BoulderSystemForSimbody	boulder_multibody(boulder, "Boulder", flap_multibody_shape);
	/** Mass properties of the consrained spot. 
	 * SimTK::MassProperties(mass, center of mass, inertia)
	 */
	SimTK::Body::Rigid boulder_info(*boulder_multibody.body_part_mass_properties_);
	/** Adding the simbody contact to bolder. */
	addBoulderContactForSimbody(boulder_info);

	SimTK::MobilizedBody::Planar boulder_body(matter.Ground(), 
		SimTK::Transform(Vec3d(B_x - BL/2.0, B_y + BH/2.0, 0.0)), 
		boulder_info, SimTK::Transform());
	cout << "SUCCESS!\n";
	/** Time steping method for multibody system.*/
	cout << "Realizing topology... ";
	SimTK::State state = MBsystem.realizeTopology();
	cout << "SUCCESS!\n";
	cout << "Selecting time stepping method... ";
	SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
	// SimTK::CPodesIntegrator integ(MBsystem, SimTK::CPodes::BDF, SimTK::CPodes::Newton);
	integ.setAccuracy(1e-3);
	integ.setAllowInterpolation(false);
	cout << "SUCCESS!\n";
	cout << "Initializing simbody state... ";
	integ.initialize(state);
	cout << "SUCCESS!\n";

	/**
	* Coupling between SimBody and SPH.
	*/
	cout << "Simbody-SPH coupling starting... ";
	solid_dynamics::TotalForceOnSolidBodyPartForSimBody
		force_on_boulder(boulder, boulder_multibody, MBsystem, boulder_body, force_on_bodies, integ);
	solid_dynamics::ConstrainSolidBodyPartBySimBody
		constraint_boulder(boulder, boulder_multibody, MBsystem, boulder_body, force_on_bodies, integ);
	cout << "SUCCESS!\n";
	/** Output. */
	cout << "Output setup... ";
	In_Output in_output(system);
	BodyStatesRecordingToLegacyVtk write_real_body_states(in_output, system.real_bodies_);
	BodyReducedQuantityRecording<solid_dynamics::TotalForceOnSolid> write_total_force_on_boulder(in_output, boulder);
	cout << "SUCCESS!\n";
	/**
	 * @brief Prepare quantities will be used once only and initial condition.
	 */
	/** offset particle position */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromBodyShape();
	boulder_particles.initializeNormalDirectionFromBodyShape();
	boulder_corrected_configuration.parallel_exec();

	write_real_body_states.writeToFile(0);
	write_total_force_on_boulder.writeToFile(0);
	/** Simulation start here. */
	/** starting time zero. */
	system.restart_step_ = 0;
	GlobalStaticVariables::physical_time_ = 0.0;
	int number_of_iterations = 0;
	int screen_output_interval = 500;
	Real End_Time = 3.0;					/**< End time. */
	Real D_Time = End_Time / 1200.0;		/**< Time stamps for output of body states. */
	Real Dt = 0.0;							/**< Default advection time step sizes. */
	Real dt = 0.0; 							/**< Default acoustic time step sizes. */
	Real total_time = 0.0;
	Real relax_time = 0.5;
	/** statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	/** Main Loop. */
	while (GlobalStaticVariables::physical_time_ < End_Time){
		Real integral_time = 0.0;
		while (integral_time < D_Time) {
			/** Acceleration due to viscous force and gravity. */
			initialize_gravity_to_fluid.parallel_exec();
			
			Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_density_by_summation.parallel_exec();
			viscous_acceleration.parallel_exec();
			/** Viscous force exerting on boulder. */
			fluid_viscous_force_on_boulder.parallel_exec(Dt);
			boulder_update_normal.parallel_exec();
			
			/** Dynamics including pressure relaxation. */
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) {
				pressure_relaxation.parallel_exec(dt);
				fluid_pressure_force_on_boulder.parallel_exec(dt);
				density_relaxation.parallel_exec(dt);
				/** solid dynamics. */
				average_velocity_and_acceleration.initialize_displacement_.parallel_exec();

				if (total_time > relax_time)
				{
					SimTK::State& state_for_update = integ.updAdvancedState();
					// Real angle = boulder_body
					// 				.getBodyRotation(state_for_update)
					// 				.convertOneAxisRotationToOneAngle(SimTK::ZAxis);
				
					SimTK::SpatialVec current_force = force_on_boulder.parallel_exec();
					force_on_bodies.clearAllBodyForces(state_for_update);
					force_on_bodies.setOneBodyForce(state_for_update, boulder_body, 
													current_force);
					integ.stepBy(dt);
					constraint_boulder.parallel_exec();
				}
				
				average_velocity_and_acceleration.update_averages_.parallel_exec(dt);

				dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integral_time += dt;
				total_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			
			if (number_of_iterations % screen_output_interval == 0) {
				cout << fixed << setprecision(9) << "N=" << number_of_iterations 
					 << "	Total Time = " << total_time 
					 << "	Physical Time = " << GlobalStaticVariables::physical_time_ 
					 << "	Dt = " << Dt << "	dt = " << dt << "\n";
			}
			number_of_iterations++;
			water_block.updateCellLinkedList();
			wall_boundary.updateCellLinkedList();
			boulder.updateCellLinkedList();
			water_block_complex.updateConfiguration();
			boulder_fluid_contact.updateConfiguration();
			write_total_force_on_boulder.writeToFile(GlobalStaticVariables::physical_time_);
		}

		tick_count t2 = tick_count::now();
		write_real_body_states.writeToFile(GlobalStaticVariables::physical_time_ / (D_Time * 1e4));
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;

	cout << fixed << setprecision(9) << "N=" << number_of_iterations 
					 << "	Total Time = " << total_time 
					 << "	Physical Time = " << GlobalStaticVariables::physical_time_ 
					 << "	Dt = " << Dt << "	dt = " << dt << "\n";
	cout << "Total wall time for computation: " << tt.seconds() << " seconds.\n";
	
	return 0;
}
