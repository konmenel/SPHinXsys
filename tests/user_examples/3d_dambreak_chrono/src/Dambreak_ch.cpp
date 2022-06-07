/* ---------------------------------------------------------------------------*
*                       SPHinXsys: 3D dambreak example                        *
* ----------------------------------------------------------------------------*
* This is the one of the basic test cases for efficient and accurate time     *
* integration scheme investigation 							  				  *
* ---------------------------------------------------------------------------*/
#include "Chrono.h"


int main()
{
	SPHSystem system(system_domain_bounds, resolution_ref);
	LogOutput fcout("Run.out");
	GlobalStaticVariables::physical_time_ = 0.0;

	WallBoundary wall_boundary(system, "Wall");
	SolidParticles wall_particles(wall_boundary);

	Box box(system, "Box");
	ElasticSolidParticles boulder_particles(box, makeShared<LinearElasticSolid>(rho0_s, poisson, Youngs_modulus));

	fcout << "Setting up Chrono..." << std::flush;
	SystemCh ch_system;
	
	std::shared_ptr<ChBody> box_ch = addBoxCh(ch_system);
	addWallsCh(ch_system);

	ch_system.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);
	ch_system.SetSolverType(ChSolver::Type::PSOR);
	ch_system.SetSolverMaxIterations(50);
	ch_system.Set_G_acc(ChVector<>(0.0, 0.0, -9.81));

	BoulderSystemForChrono box_for_chrono(box, "Box_chrono", box.body_shape_);
	ConstrainSolidBodyPartByChrono box_constain(box, box_for_chrono, box_ch);

	fcout << "OK!" << endl;

	// Output system
	In_Output in_output(system);
	BodyStatesRecordingToLegacyVtk write_body_states(in_output, system.real_bodies_);

	size_t step = 0;
	Real dt = 0.001;
	Real end_time = 10.0;
	Real real_time = 0.0;
	Real out_time = 0.02;
	size_t report_steps = 500;

	write_body_states.writeToFile(0);

	while (real_time < end_time) {
		Real integration_time = 0.0;
		while (integration_time < out_time) {
			if (step % report_steps == 0) {
				fcout << "Step = " << step << "\treal time = " << real_time << endl;
			}

			ch_system.DoStepDynamics(dt);
			box_constain.parallel_exec();

			integration_time += dt;
			real_time += dt;
			step++;
		}
		write_body_states.writeToFile(step);
	}


	return 0;
}
