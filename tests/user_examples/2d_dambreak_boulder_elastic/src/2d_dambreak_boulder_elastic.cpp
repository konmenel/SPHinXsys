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
	/** Boulder system. Body, material and particle container. */
	Boulder boulder(system, "Boulder");
	ElasticSolidParticles boulder_particles(boulder, makeShared<BoulderMaterial>());
	/** topology */
	BodyRelationInner boulder_inner(boulder);
	ComplexBodyRelation water_block_complex(water_block, { &wall_boundary, &boulder });
	BodyRelationContact boulder_fluid_contact(boulder, { &water_block });
	SolidBodyRelationContact boulder_wall_contact(boulder, { &wall_boundary });
	/**
	 * Methods only used only once
	 */
	/** corrected strong configuration. */
	solid_dynamics::CorrectConfiguration boulder_corrected_configuration(boulder_inner);
	/** 
	 * Methods used for time stepping
	 */
	/** Time step initialization, add gravity. */
	TimeStepInitialization 	initialize_gravity_to_fluid(water_block, gravity);
	TimeStepInitialization 	initialize_gravity_to_boulder(boulder, gravity);
	/** Evaluation of density by summation approach. */
	fluid_dynamics::DensitySummationFreeSurfaceComplex update_density_by_summation(water_block_complex);
	/** time step size without considering sound wave speed. */
	fluid_dynamics::AdvectionTimeStepSize	get_fluid_advection_time_step_size(water_block, U_f);
	/** time step size with considering sound wave speed. */
	fluid_dynamics::AcousticTimeStepSize	get_fluid_time_step_size(water_block);
	solid_dynamics::AcousticTimeStepSize 	boulder_get_time_step_size(boulder);
	/** pressure relaxation using verlet time stepping. */
	solid_dynamics::StressRelaxationFirstHalf boulder_stress_relaxation_first_half(boulder_inner);
	solid_dynamics::StressRelaxationSecondHalf boulder_stress_relaxation_second_half(boulder_inner);
	fluid_dynamics::PressureRelaxationRiemannWithWall	pressure_relaxation(water_block_complex);
	fluid_dynamics::DensityRelaxationRiemannWithWall	density_relaxation(water_block_complex);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAccelerationWithWall viscous_acceleration(water_block_complex);
	fluid_dynamics::TransportVelocityCorrectionComplex	transport_velocity_correction(water_block_complex);
	CombinedInteractionDynamics viscous_acceleration_and_transport_correction(viscous_acceleration, transport_velocity_correction);
	/** Fluid force on boulder. */
	solid_dynamics::FluidPressureForceOnSolid fluid_pressure_force_on_boulder(boulder_fluid_contact);
	/** Fluid viscous force on boulder. */
	solid_dynamics::FluidViscousForceOnSolid fluid_viscous_force_on_boulder(boulder_fluid_contact);
	/** Contact force on boulder. */
	solid_dynamics::ContactDensitySummation boulder_update_contact_density(boulder_wall_contact);
	// solid_dynamics::ContactForce contact_force_on_boulder(boulder_wall_contact);
	solid_dynamics::ContactForceWithWall contact_force_on_boulder(boulder_wall_contact);
	/** Damping */
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec2d>>
		boulder_damping(boulder_inner, 0.5, "Velocity", physical_viscosity);
	/** average velocity for boulder. */
	solid_dynamics::AverageVelocityAndAcceleration	average_velocity_and_acceleration(boulder);
	solid_dynamics::UpdateElasticNormalDirection 	boulder_update_normal(boulder);
	
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
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	wall_particles.initializeNormalDirectionFromBodyShape();
	boulder_particles.initializeNormalDirectionFromBodyShape();
	boulder_corrected_configuration.parallel_exec();

	write_real_body_states.writeToFile(GlobalStaticVariables::physical_time_);
	write_total_force_on_boulder.writeToFile(GlobalStaticVariables::physical_time_);
	/** Simulation start here. */
	/** starting time zero. */
	system.restart_step_ = 0;
	GlobalStaticVariables::physical_time_ = 0.0;
	int number_of_iterations = 0;
	int screen_output_interval = 500;
	Real End_Time = 2.0;					/**< End time. */
	Real D_Time = End_Time / 1200.0;		/**< Time stamps for output of body states. */
	Real Dt = 0.0;							/**< Default advection time step sizes. */
	Real dt = 0.0; 							/**< Default acoustic time step sizes. */
	Real dt_s = 0.0;						/**< Default acoustic time step sizes for solid. */
	size_t inner_ite_dt = 0;
	size_t inner_ite_dt_s = 0;
	Real total_time = 0.0;
	/** statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/** Main Loop. */
	write_real_body_states.writeToFile(0);
	while (GlobalStaticVariables::physical_time_ < End_Time){
		Real integral_time = 0.0;
		
		while (integral_time < D_Time) {
			/** Acceleration due to viscous force and gravity. */
			initialize_gravity_to_fluid.parallel_exec();
			initialize_gravity_to_boulder.parallel_exec();
			
			Dt = get_fluid_advection_time_step_size.parallel_exec();
			update_density_by_summation.parallel_exec();
			viscous_acceleration_and_transport_correction.parallel_exec();
			
			fluid_viscous_force_on_boulder.parallel_exec();
			boulder_update_normal.parallel_exec();
			
			/** Dynamics including pressure relaxation. */
			Real relaxation_time = 0.0;
			inner_ite_dt = 0;
			while (relaxation_time < Dt) {
				dt = SMIN(get_fluid_time_step_size.parallel_exec(), Dt);

				pressure_relaxation.parallel_exec(dt);
				boulder_update_contact_density.parallel_exec();
				contact_force_on_boulder.parallel_exec();
				fluid_pressure_force_on_boulder.parallel_exec();
				boulder_damping.parallel_exec(dt);
				density_relaxation.parallel_exec(dt);
				/** solid dynamics. */
				inner_ite_dt_s = 0;
				Real dt_s_sum = 0.0;
				average_velocity_and_acceleration.initialize_displacement_.parallel_exec();
				while (dt_s_sum < dt) {
					dt_s = SMIN(boulder_get_time_step_size.parallel_exec(), dt - dt_s_sum);
					boulder_stress_relaxation_first_half.parallel_exec(dt_s);
					boulder_stress_relaxation_second_half.parallel_exec(dt_s);
					dt_s_sum += dt_s;
					inner_ite_dt_s++;
				}
				average_velocity_and_acceleration.update_averages_.parallel_exec(dt);

				// dt = get_fluid_time_step_size.parallel_exec();
				relaxation_time += dt;
				integral_time += dt;
				total_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				inner_ite_dt++;
			}
			
			if (number_of_iterations % screen_output_interval == 0) {
				cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations 
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
			boulder_wall_contact.updateConfiguration();
			write_total_force_on_boulder.writeToFile(number_of_iterations);
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
