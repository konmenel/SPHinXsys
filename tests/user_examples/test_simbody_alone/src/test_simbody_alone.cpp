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
	/** Multi-body system. */
	cout << "Simbody setup... ";
	/** set up the multi body system. */
	MultibodySystem system;
	/** the bodies or matter of the system. */
	SimbodyMatterSubsystem 	matter(system);
	/** the forces of the system. */
	GeneralForceSubsystem 	forces(system);
	/** the contact system. */
	ContactTrackerSubsystem  tracker(system);
	CompliantContactSubsystem contactForces(system, tracker);
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
	Force::UniformGravity sim_gravity(forces, matter, Vec3(0.0, -gravity_g, 0.0));
	/** discrete forces acting on the bodies. */
	// Force::DiscreteForces force_on_bodies(forces, matter);
	cout << "SUCCESS!\n";
	
	/** Wall contacts */
	cout << "Adding wall contacts... ";
	/** contanct clique for wall */
	ContactCliqueId clique = ContactSurface::createNewContactClique();
	
	addSimbodyWallContacts(matter, clique);
	addWallContactsDecor(matter);

	// Cliff body
	MassProperties cliff_properties(
		1e7, Vec3(0.0), UnitInertia::brick(cliff_hlx, cliff_hly, 0.0));
	Body::Rigid cliff_body(cliff_properties);
	// addCliffContactForSimbody(cliff_body, clique);
	addCliffContactForSimbody(matter, clique);
	addCliffContactDecor(matter);

	// MobilizedBody::Translation cliff(matter.Ground(), Transform(cliff_center),
	// 		cliff_body, Transform());
	// Constraint::ConstantSpeed(cliff, MobilizerUIndex(0), 0.0);
	// Constraint::ConstantSpeed(cliff, MobilizerUIndex(1), 0.0);
	// Constraint::ConstantSpeed(cliff, MobilizerUIndex(2), 0.0);
	cout << "SUCCESS!\n";

	const int nsurf = matter.Ground().getBody().getNumContactSurfaces();
	cout << "Number of Contact Surfaces added to ground = " << nsurf << "\n";
	
	/** mass properties of boulder. */
	cout << "Adding boulder... ";
	MassProperties boulder_properties(
				boulder_mass, 
				SimTK::Vec3(0.0), 
				SimTK::UnitInertia::brick(BL/2.0, BH/2.0, 0.0));
	/** Mass properties of the consrained spot. 
	 * MassProperties(mass, center of mass, inertia)
	 */
	Body::Rigid boulder_info(boulder_properties);
	/** Adding the simbody contact to bolder. */
	// ContactGeometry::Brick boulder_sim_geometry(Vec3(BL/2.0, BH/2.0, 0.01));
	// boulder_info.addContactSurface(Transform(),
    //     ContactSurface(boulder_sim_geometry,
    //                    ContactMaterial(fK, fDis, fFac, fFac, fVis)));
	addBoulderContactForSimbody(boulder_info);
	addBoulderDecor(boulder_info);

	MobilizedBody::Free boulder_body(matter.Ground(), 
		Transform(Vec3(B_x - BL/2.0, B_y + BH/2.0, 0.0)), 
		boulder_info, Transform());

	// boulder 2
	// Body::Rigid boulder_info2(boulder_properties);
	// addBoulderContactForSimbody(boulder_info2);
	// addBoulderDecor(boulder_info2, Red);

	// MobilizedBody::Free boulder_body2(matter.Ground(), 
	// 	Transform(Vec3(B_x - BL/2.0, B_y - 0.05 + BH/2.0, 0.0)), 
	// 	boulder_info2, Transform());
	cout << "SUCCESS!\n";

	/** Time steping method for multibody system.*/
	cout << "Initializing simbody visualizer... ";
	Visualizer viz(system);
    viz.setMode(Visualizer::RealTime);
    viz.setDesiredBufferLengthInSec(0.1);
    viz.setDesiredFrameRate(FPS);
    viz.setGroundHeight(-BW);
    viz.setShowShadows(false);
	viz.setBackgroundType(Visualizer::SolidColor);
	viz.setBackgroundColor(White);
	cout << "SUCCESS!\n";

	system.addEventReporter(new MyReporter(system, contactForces, ReportInterval));
    system.addEventReporter(new Visualizer::Reporter(viz, ReportInterval));

	cout << "Realizing topology... ";
	system.realizeTopology();
	cout << "SUCCESS!\n";
    
    State state = system.getDefaultState();
	viz.report(state);

	cout << "Selecting time stepping method... ";
	// RungeKuttaMersonIntegrator integ(system);
	CPodesIntegrator integ(system, CPodes::BDF, CPodes::Newton);
	integ.setAccuracy(1e-4);
	cout << "SUCCESS!\n";

	cout << "Initializing simbody state... ";
	TimeStepper ts(system, integ);
	ts.initialize(state);
	cout << "SUCCESS!\n";

	/** Simulation start here. */
	/** starting time zero. */
	GlobalStaticVariables::physical_time_ = 0.0;
	Real End_Time = 1.0;					/**< End time. */
	// Real Dt = 0.001;

	std::string user_in;
	cout << "Press `Enter` to start...";
	getchar();

	// while (GlobalStaticVariables::physical_time_ < End_Time) {
	// 	integ.stepBy(Dt);
	// 	GlobalStaticVariables::physical_time_ += Dt;
	// 	cout << GlobalStaticVariables::physical_time_ << "\n";
	// }
	ts.stepTo(End_Time);
	cout << "Done!\n";
    viz.dumpStats(cout);

    viz.setMode(Visualizer::PassThrough);
	double speed = 0.1;
	while (true) {
		cout << "Again? ";
		getchar();
		for (double i=0; i < (int)saveEm.size(); i += speed ) {
            viz.report(saveEm[(int)i]);
        }
	}
	return 0;
}
