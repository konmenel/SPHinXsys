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
	/** The wall boundary, body and particles container. */
	WallBoundary wall_boundary(system, "Wall");
	SolidParticles wall_particles(wall_boundary);
	/** Boulder system. Body, material and particle container. */
	Boulder boulder(system, "Boulder");
	ElasticSolidParticles boulder_particles(boulder, makeShared<BoulderMaterial>());
	/** topology */
	BodyRelationInner boulder_inner(boulder);
	SolidBodyRelationContact boulder_wall_contact(boulder, { &wall_boundary });
	/**
	 * Methods only used only once
	 */
	/** corrected strong configuration. */
	solid_dynamics::CorrectConfiguration boulder_corrected_configuration(boulder_inner);
	/** 
	 * Methods used for time stepping
	 */
	/** Contact force on boulder. */
	solid_dynamics::ContactDensitySummation boulder_update_contact_density(boulder_wall_contact);
	// solid_dynamics::ContactForce contact_force_on_boulder(boulder_wall_contact);
	solid_dynamics::ContactForceWithWall contact_force_on_boulder(boulder_wall_contact);
	/** Damping*/
	// DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec2d>>
	// 	boulder_damping(boulder_inner, 0.5, "Velocity", physical_viscosity);
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
	/** mass properties of boulder. */
	cout << "SUCCESS!\n";

	/** Add gravity to mb body. */
	cout << "Adding forces... ";
	SimTK::Force::UniformGravity sim_gravity(forces, matter, SimTK::Vec3(0.0, -gravity_g, 0.0));
	/** discrete forces acting on the bodies. */
	SimTK::Force::DiscreteForces force_on_bodies(forces, matter);
	cout << "SUCCESS!\n";

	cout << "Adding boulder... ";
	MultiPolygonShape boulder_shape(CreateBoulderMultiShape());
	BoulderSystemForSimbody boulder_multibody(boulder, "Boulder", boulder_shape);
	/** Mass properties of the consrained spot. 
	 * SimTK::MassProperties(mass, center of mass, inertia)
	 */
	SimTK::Body::Rigid      boulder_info(*boulder_multibody.body_part_mass_properties_);
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
	SimTK::MobilizedBody::Planar boulder_body(matter.Ground(),
		SimTK::Transform(Vec3d(B_x - BL/2.0, B_y + BH/2.0, 0.0)), 
		boulder_info, SimTK::Transform(Vec3d(0.0, 0.0, 0.0)));
	cout << "SUCCESS!\n";

	/** Time steping method for multibody system.*/
	cout << "Initializing simbody... ";
	SimTK::State state = MBsystem.realizeTopology();
	SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
	integ.setAccuracy(1e-3);
	integ.setAllowInterpolation(false);
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
	BodyStatesRecordingToVtp 		write_real_body_states(in_output, system.real_bodies_);
	RegressionTestDynamicTimeWarping<BodyReducedQuantityRecording<solid_dynamics::TotalForceOnSolid>> 
		write_total_force_on_boulder(in_output, boulder);
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
	// BoulderInitialCondition boulder_initial_vel(boulder);
	// boulder_initial_vel.parallel_exec();
	boulder_body.setOneU(state, 2, 1.0);

	write_real_body_states.writeToFile(0);
	write_total_force_on_boulder.writeToFile(0);
	/** Simulation start here. */
	/** starting time zero. */
	system.restart_step_ = 0;
	GlobalStaticVariables::physical_time_ = 0.0;
	int number_of_iterations = 0;
	int screen_output_interval = 500;
	Real End_Time = 1.0;					/**< End time. */
	Real D_Time = End_Time / 1500.0;		/**< Time stamps for output of body states. */
	Real Dt = 1e-4;							/**< Timestep */
	Real total_time = 0.0;
	/** statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/** Main Loop. */
	write_real_body_states.writeToFile(GlobalStaticVariables::physical_time_);
	while (GlobalStaticVariables::physical_time_ < End_Time) {
		Real integral_time = 0.0;
		while (integral_time < D_Time) {
			boulder_update_normal.parallel_exec();
			boulder_update_contact_density.parallel_exec();
			contact_force_on_boulder.parallel_exec();
			// boulder_damping.parallel_exec(Dt);
			/** solid dynamics. */
			average_velocity_and_acceleration.initialize_displacement_.parallel_exec();

			SimTK::State& state_for_update = integ.updAdvancedState();
			force_on_bodies.clearAllBodyForces(state_for_update);
			force_on_bodies.setOneBodyForce(state_for_update, boulder_body, 
											force_on_boulder.parallel_exec());
			integ.stepBy(Dt);
			constraint_boulder.parallel_exec();
			
			average_velocity_and_acceleration.update_averages_.parallel_exec(Dt);

			integral_time += Dt;
			total_time += Dt;
			GlobalStaticVariables::physical_time_ += Dt;
			
				
			if (number_of_iterations % screen_output_interval == 0) {
				cout << fixed << setprecision(9) << "N=" << number_of_iterations 
						<< "	Total Time = " << total_time 
						<< "	Physical Time = " << GlobalStaticVariables::physical_time_ 
						<< "	Dt = " << Dt << "\n";
			}
			number_of_iterations++;
			wall_boundary.updateCellLinkedList();
			boulder.updateCellLinkedList();
			boulder_wall_contact.updateConfiguration();
			write_total_force_on_boulder.writeToFile(GlobalStaticVariables::physical_time_);
		}
		tick_count t2 = tick_count::now();
		write_real_body_states.writeToFile(number_of_iterations);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds.\n";
	
	return 0;
}
