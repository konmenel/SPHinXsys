/* ---------------------------------------------------------------------------*
*                       SPHinXsys: 3D dambreak example                        *
* ----------------------------------------------------------------------------*
* This is the one of the basic test cases for efficient and accurate time     *
* integration scheme investigation 							  				  *
* ---------------------------------------------------------------------------*/
#include "test_vtk.h"


//the main program
int main()
{
	cout << "System Endianness: " << ((Endian::getSystemEndianness()==Endianness::little) ? "little" : "big") << endl;
	double position[9];
	int ids[3] = { 0, 1, 2 };
	int unsorted_ids[3] = { 1, 0, 2 };
	double vel[9] = { 0.0, 0.0, 0.0,
					  1.0, 1.0, 1.0,
					  2.0, 2.0, 2.0 };
	double acc[9] = { 1.0, 1.0, 1.0,
					  0.0, 0.0, 0.0,
					  2.0, 2.0, 2.0 };

	// initialize the position with random values
	srand(42);
	for (int i = 0; i < 9; i++)
		position[i] = rand() / (double)RAND_MAX;
	//create a new VTK file
	ofstream vtk_file("test_vtk.vtk");
	//write the header
	vtk_file << "# vtk DataFile Version 3.0" << endl;
	vtk_file << "Test data" << endl;
	vtk_file << "BINARY" << endl;
	vtk_file << "DATASET POLYDATA" << endl;
	vtk_file << "POINTS 3 float" << endl;
	//write the data
	for (int i = 0; i < 9; i++) {
		float position_float = static_cast<float>(position[i]);
		Endian::writeDataReverseEndianness(vtk_file, &position_float, sizeof(float), sizeof(float));
	}
	vtk_file << endl;
	vtk_file << "POINT_DATA 3" << endl;
	vtk_file << "SCALARS Particle_ID unsigned_int" << endl;
	vtk_file << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i < 3; i++) {
		unsigned int id = ids[i];
		Endian::writeDataReverseEndianness(vtk_file, &id, sizeof(unsigned int), sizeof(unsigned int));
	}
	vtk_file << endl;
	vtk_file << "FIELD FieldData 2" << endl;
	vtk_file << "Vel 3 3 float" << endl;
	for (int i = 0; i < 9; i++) {
		float vel_float = static_cast<float>(vel[i]);
		Endian::writeDataReverseEndianness(vtk_file, &vel_float, sizeof(float), sizeof(float));
	}
	vtk_file << endl;
	vtk_file << "Acc 3 3 float" << endl;
	for (int i = 0; i < 9; i++) {
		float vel_float = static_cast<float>(acc[i]);
		Endian::writeDataReverseEndianness(vtk_file, &vel_float, sizeof(float), sizeof(float));
	}
	vtk_file << endl;
	vtk_file.close();

	// read the data
	ifstream vtk_file_read("test_vtk.vtk", ios::binary);
	// skip 5 lines (the header)
	string header;
	for (int i = 0; i < 5; i++) {
		getline(vtk_file_read, header);
	}
	// read the data
	float position_read[9];
	for (int i = 0; i < 9; i++) {
		for (int j = sizeof(float)-1; j >= 0; j--) {
			vtk_file_read.read(((char *)&position_read[i])+j, 1);
		}
	}
	vtk_file_read.close();
	// check the data
	cout << "Checking the data..." << endl;
	for (int i = 0; i < 9; i++) {
		cout << position_read[i] << " " << position[i] << endl;
	}

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

	/** topology */
	ComplexBodyRelation water_block_complex(water_block, {&wall_boundary});

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
	fluid_dynamics::PressureRelaxationRiemannWithWall pressure_relaxation(water_block_complex);
	fluid_dynamics::DensityRelaxationRiemannWithWall density_relaxation(water_block_complex);
	/** Computing viscous acceleration. */
	fluid_dynamics::ViscousAccelerationWithWall viscous_acceleration(water_block_complex);
	/** Impose transport velocity. */
	fluid_dynamics::TransportVelocityCorrectionComplex transport_velocity_correction(water_block_complex);
	/** viscous acceleration and transport velocity correction can be combined because they are independent dynamics. */
	CombinedInteractionDynamics viscous_acceleration_and_transport_correction(viscous_acceleration, transport_velocity_correction);

	//-----------------------------------------------------------------------------
	//outputs
	//-----------------------------------------------------------------------------
	In_Output in_output(system);
	BodyStatesRecordingToLegacyVtk write_water_block_states(in_output, system.real_bodies_);
	BodyStatesRecordingToVtp write_water_block_states_vtp(in_output, system.real_bodies_);
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
	water_block.setNewlyUpdated();
	write_water_block_states_vtp.writeToFile(0);

	size_t number_of_iterations = system.restart_step_;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 1.6;
	//time step size for output file
	Real D_Time = 0.1;
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
			// transport_velocity_correction.parallel_exec();
			viscous_acceleration_and_transport_correction.parallel_exec();

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
		}


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
