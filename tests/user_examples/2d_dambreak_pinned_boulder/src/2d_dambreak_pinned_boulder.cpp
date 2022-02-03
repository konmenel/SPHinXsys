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

// TODO: Check simbody constrains
// TODO: change all flap related to boulder related
// TODO: add dumping
// TODO: add observers

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
	BoulderMaterial* boulder_material = new BoulderMaterial();
	ElasticSolidParticles boulder_particles(boulder, boulder_material);
	/** Pressure probe on Flap. */
	// FlapObserver* observer = new FlapObserver(system, "FlapObserver");
	// BaseParticles 	observer_particles(observer);
	/** topology */
	// InnerBodyRelation* water_block_inner 	= new InnerBodyRelation(water_block);
	InnerBodyRelation* boulder_inner 			= new InnerBodyRelation(boulder);
	ComplexBodyRelation* water_block_complex = new ComplexBodyRelation(water_block, { wall_boundary, boulder });
	ContactBodyRelation* boulder_fluid_contact 		= new ContactBodyRelation(boulder, { water_block });
	SolidContactBodyRelation* boulder_wall_contact 		= new SolidContactBodyRelation(boulder, { wall_boundary });
	// ContactBodyRelation* observer_contact_with_water = new ContactBodyRelation(observer, { water_block });
	// ContactBodyRelation* observer_contact_with_flap  = new ContactBodyRelation(observer, {flap});
	/**
	 * Methods only used only once
	 */
	/** corrected strong configuration. */
	solid_dynamics::CorrectConfiguration boulder_corrected_configuration(boulder_inner);
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
	/** Inflow boundary condition. */
	// fluid_dynamics::DampingBoundaryCondition	damping_wave(water_block, new DampingBuffer(water_block, "DampingBuffer"));
	/** Fluid force on boulder. */
	solid_dynamics::FluidPressureForceOnSolid fluid_pressure_force_on_boulder(boulder_fluid_contact);
	/** Fluid viscous force on boulder. */
	solid_dynamics::FluidViscousForceOnSolid fluid_viscous_force_on_boulder(boulder_fluid_contact);
	/** Contact force on boulder. */
	solid_dynamics::ContactForce contact_force_on_boulder(boulder_wall_contact);
	// solid_dynamics::ContactForceFromFriction friction_force_on_boulder(boulder_wall_contact);
	/** average velocity for boulder. */
	solid_dynamics::AverageVelocityAndAcceleration	average_velocity_and_acceleration(boulder);
	solid_dynamics::UpdateElasticNormalDirection 	boulder_update_normal(boulder);
	/** constrain region of the part of wall boundary. */
	// WaveMaking wave_making(wall_boundary, new WaveMaker(wall_boundary, "WaveMaker"));
	
	/** Multi-body system. */
	std::cout << "Simbody setup... ";
	/** set up the multi body system. */
	SimTK::MultibodySystem MBsystem;
	/** the bodies or matter of the system. */
	SimTK::SimbodyMatterSubsystem 	matter(MBsystem);
	/** the forces of the system. */
	SimTK::GeneralForceSubsystem 	forces(MBsystem);
	/** mass properties of boulder. */
	std::cout << "SUCCESS!\n";
	std::cout << "Adding boulder... ";
	BoulderSystemForSimbody*  	boulder_multibody = new BoulderSystemForSimbody(boulder, "Boulder", rho0_s);
	/** Mass properties of the consrained spot. 
	 * SimTK::MassProperties(mass, center of mass, inertia)
	 */
	SimTK::Body::Rigid      boulder_info(*boulder_multibody->body_part_mass_properties_);
	/** 
	 * @brief   Pin (MobilizedBody &parent, const Transform &X_PF, const Body &bodyInfo, const 
	 					Transform &X_BM, Direction=Forward)1
	 * @details Create a Pin mobilizer between an existing parent (inboard) body P and 
	 * 			a new child (outboard) body B created by copying the given bodyInfo into 
	 *			a privately-owned Body within the constructed MobilizedBody object.
	 * @param[in] inboard(SimTK::Vec3) Defines the location of the joint point relative to the parent body.
	 * @param[in] outboard(SimTK::Vec3) Defines the body's origin location to the joint point. 
	 * @note	The body's origin location can be the mass center, the the center of mass should be Vec3(0)
	 * 			in SimTK::MassProperties(mass, com, inertia)
	 */
	SimTK::MobilizedBody::Pin boulder_body(matter.Ground(), SimTK::Transform(SimTK::Vec3(B_x - BW/2 , B_y + BH, 0.0)), 
		boulder_info, SimTK::Transform(SimTK::Vec3(0.0, 0.0, 0.0)));
	// SimTK::MobilizedBody::Free boulder_body(matter.Ground(), SimTK::Transform(SimTK::Vec3(B_x, 0.0, 0.0)), 
	// 	boulder_info, SimTK::Transform(SimTK::Vec3(0.0, 0.0, 0.0)));
	std::cout << "SUCCESS!\n";
	/** set the default angle of the pin. */
	// boulder_body.setDefaultAngle(0);
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
	std::cout << "Adding forces... ";
	SimTK::Force::UniformGravity sim_gravity(forces, matter, SimTK::Vec3(0.0, -9.81, 0.0));
	/** discrete forces acting on the bodies. */
	SimTK::Force::DiscreteForces force_on_bodies(forces, matter);
	std::cout << "SUCCESS!\n";
	/**
	 * Add a linear damping force to the mobilized body.
	 * @class SimTK::Force::MobilityLinearDamper::MobilityLinearDamper( 	
	 * @param[in]	GeneralForceSubsystem &  	forces,
	 * @param[in]	const MobilizedBody &  	mobod,
     * @param[in]	MobilizerUIndex  whichU, e.g., MobilizerUIndex(0)
	 * @param[in]	Real  	Dampingconstant ) 
	 * Here, The damping constant c is provided, with the generated force being -c*u where u is the mobility's generalized speed.	
	 */
	// SimTK::Force::MobilityLinearDamper linear_damper(forces, pin_spot, SimTK::MobilizerUIndex(0), 20.0);
	/** Time steping method for multibody system.*/
	std::cout << "Realizing topology... ";
	SimTK::State state = MBsystem.realizeTopology();
	std::cout << "SUCCESS!\n";
	std::cout << "Selecting time stepping method... ";
	SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
	integ.setAccuracy(1e-3);
	integ.setAllowInterpolation(false);
	std::cout << "SUCCESS!\n";
	std::cout << "Initializing simbody state... ";
	integ.initialize(state);
	std::cout << "SUCCESS!\n";
	/**
	* Coupling between SimBody and SPH.
	*/
	std::cout << "Simbody-SPH coupling starting... ";
	solid_dynamics::TotalForceOnSolidBodyPartForSimBody
		force_on_boulder(boulder, boulder_multibody, MBsystem, boulder_body, force_on_bodies, integ);
	solid_dynamics::ConstrainSolidBodyPartBySimBody
		constraint_boulder(boulder, boulder_multibody, MBsystem, boulder_body, force_on_bodies, integ);
	std::cout << "SUCCESS!\n";
	/** Output. */
	std::cout << "Output setup... ";
	In_Output in_output(system);
	WriteBodyStatesToVtu 		write_real_body_states(in_output, system.real_bodies_);
	WriteTotalForceOnSolid      write_total_force_on_boulder(in_output, boulder);
	WriteSimBodyPinData			write_boulder_pin_data(in_output, integ, boulder_body);
	std::cout << "SUCCESS!\n";
	/**
	 * @brief Prepare quantities will be used once only and initial condition.
	 */
	/** offset particle position */
	// flap_particles.offsetInitialParticlePosition(offset);
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromGeometry();
	boulder_particles.initializeNormalDirectionFromGeometry();
	boulder_corrected_configuration.parallel_exec();

	write_real_body_states.WriteToFile(GlobalStaticVariables::physical_time_);
	write_total_force_on_boulder.WriteToFile(GlobalStaticVariables::physical_time_);
	write_boulder_pin_data.WriteToFile(GlobalStaticVariables::physical_time_);
	/** Simulation start here. */
	/** starting time zero. */
	system.restart_step_ = 0;
	GlobalStaticVariables::physical_time_ = 0.0;
	int number_of_iterations = 0;
	int screen_output_interval = 500;
	Real End_Time = 2.0;					/**< End time. */
	Real D_Time = End_Time / 1000.0;		/**< Time stamps for output of body states. */
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
			/** Viscous force exerting on boulder. */
			fluid_viscous_force_on_boulder.parallel_exec(Dt);
			contact_force_on_boulder.parallel_exec(Dt);
			boulder_update_normal.parallel_exec();
			
			/** Dynamics including pressure relaxation. */
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) {
				pressure_relaxation.parallel_exec(dt);
				fluid_pressure_force_on_boulder.parallel_exec(dt);
				density_relaxation.parallel_exec(dt);
				/** solid dynamics. */
				average_velocity_and_acceleration.initialize_displacement_.parallel_exec();

				SimTK::State& state_for_update = integ.updAdvancedState();
				// const SimTK::Rotation& boulder_rotation = boulder_body.getBodyRotation(state_for_update);
				// Real angle = boulder_rotation.convertOneAxisRotationToOneAngle(SimTK::CoordinateAxis::CoordinateAxis(2));
				
				Real angle = boulder_body.getAngle(state_for_update); 	// for SimTK::MobilizedBody::Pin
				force_on_bodies.clearAllBodyForces(state_for_update);
				force_on_bodies.setOneBodyForce(state_for_update, boulder_body, force_on_boulder.parallel_exec(angle));
				integ.stepBy(dt);
				constraint_boulder.parallel_exec();
				// wave_making.parallel_exec(dt);
				
				average_velocity_and_acceleration.update_averages_.parallel_exec(dt);
				// interpolation_observer_position.parallel_exec();

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
			water_block->updateCellLinkedList();
			wall_boundary->updateCellLinkedList();
			boulder->updateCellLinkedList();
			water_block_complex->updateConfiguration();
			boulder_fluid_contact->updateConfiguration();
			boulder_wall_contact->updateConfiguration();
			write_total_force_on_boulder.WriteToFile(GlobalStaticVariables::physical_time_);
			write_boulder_pin_data.WriteToFile(GlobalStaticVariables::physical_time_ );
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
