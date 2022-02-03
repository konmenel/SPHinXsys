 /**
  * @brief 	SPHinXsys Library.
  */
#pragma once
#include "sphinxsys.h"
  /**
 * @brief Namespace cite here.
 */
using namespace SimTK;
using std::cout;
using std::endl;
using SPH::GlobalStaticVariables;


Array_<State> saveEm;
const Real FPS = 200.0;
const Real ReportInterval = 1.0 / FPS;

/**
 * @brief Basic geometry parameters and numerical setup.
 */
const Real DL = 3.0; 								/**< Tank length [m]. */
const Real DH = 0.5; 								/**< Tank height [m]. */
const Real wall_H = 0.2;							/**< Verical wall height [m] */
const Real wall_position = 2.0;						/**< Position of the verical wall [m] (x direction) */ 
const Real BL = 2.0e-2;								/**< Boulder lenght [m]. */
const Real BH = 1.5e-2;								/**< Boulder height [m]. */
const Real B_x = wall_position + 0.5;				/**< Boulder initial position x-axis (right edge) [m]. */
const Real B_y = wall_H + 0.1;						/**< Boulder initial position y-axis (right edge) [m]. */
const Real particle_spacing_ref = BH / 6.0; 		/**< Initial reference particle spacing. */
const Real BW = particle_spacing_ref * 6.0; 		/**< Extending width for BCs. */

const Real gravity_g = 9.81;		/**< Gravity force of fluid. [m/s^2] */

/**
 * @brief Material properties of the solid.
 */
const Real rho0_s = 2.8e3;								/**< Boulder Density [kg/m^3] from paper. */
const Real boulder_vol = BL * BH;						/**< Boulder Volume [m^3]. (1.5 x 2.0 x 3.0 cm) */
const Real boulder_mass = rho0_s * boulder_vol;			/**< Boulder Mass [kg]. */
const Real poisson = 0.25;								/**< Poisson's ratio. */
const Real Youngs_modulus = 73e9;						/**< Young's modulus [Pa]. */
// const Real physical_viscosity = 1000000.0;

/**
 * @brief Cliff dimensions.
 */
const Real cliff_hlx = (DL - wall_position) / 2.0;
const Real cliff_hly = wall_H / 2.0;
const Vec3 cliff_center(wall_position + cliff_hlx, cliff_hly, 0.0);

/**
 * @brief Material properties of the simbody contact.
 */
const Real fK = ContactMaterial 							/**< Stiffness Coefficient [Pa] */
		::calcPlaneStrainStiffness(Youngs_modulus, poisson);
const Real fDis = 0.01; // to turn off dissipation
const Real fFac = 0.5; // to turn off friction
const Real fVis = 0.0; // to turn off viscous friction


void addSimbodyWallContacts(SimbodyMatterSubsystem& matter, 
		const ContactCliqueId& clique)
{
	const ContactGeometry::HalfSpace half_space;
	const ContactMaterial material(fK, fDis, fFac, fFac, fVis);
	// Left wall
	const Rotation R_left(Pi, ZAxis);
	matter.Ground().updBody().addContactSurface(
		Transform(R_left),
        ContactSurface(half_space, material)
                       .joinClique(clique));
	// Floor
	const Rotation R_floor(-Pi/2.0, ZAxis);
	matter.Ground().updBody().addContactSurface(
		Transform(R_floor),
        ContactSurface(half_space, material)
                       .joinClique(clique));
	// Ceiling
	const Rotation R_ceiling(Pi/2.0, ZAxis);
	matter.Ground().updBody().addContactSurface(
		Transform(R_ceiling, Vec3(0.0, DH, 0.0)),
        ContactSurface(half_space, material)
                       .joinClique(clique));
	// Right wall
	matter.Ground().updBody().addContactSurface(Vec3(DL, 0.0, 0.0),
        ContactSurface(half_space, material)
                       .joinClique(clique));
}


void addCliffContactForSimbody(Body& cliff, const ContactCliqueId& clique) 
{
	return;
	const ContactMaterial material(fK, fDis, fFac, fFac, fVis);

	// Contact surface
	PolygonalMesh brick_mesh;
	brick_mesh = PolygonalMesh::createBrickMesh(Vec3(cliff_hlx-0.001, cliff_hly+0.001, 0.1), 0);
	ContactGeometry::TriangleMesh cliff_geometry(brick_mesh);
	
	cliff.addContactSurface(Transform(),
        ContactSurface(cliff_geometry, material, 0.001));
}


void addCliffContactForSimbody(SimbodyMatterSubsystem& matter,
		const ContactCliqueId& clique) 
{
	const ContactMaterial material(fK, fDis, fFac, fFac, fVis);
	
	// Contact surface
	PolygonalMesh brick_mesh;
	brick_mesh = PolygonalMesh::createBrickMesh(Vec3(cliff_hlx-0.001, cliff_hly-0.001, 0.1), 0);
	ContactGeometry::TriangleMesh cliff_geometry(brick_mesh);

	matter.Ground().updBody().addContactSurface(Transform(cliff_center),
        ContactSurface(cliff_geometry, material, 0.001).joinClique(clique));
}


void addBoulderContactForSimbody(Body::Rigid& boulder_body)
{
	const ContactMaterial material(fK, fDis, fFac, fFac, fVis);
	
	// Contact surface
	PolygonalMesh brick_mesh;
	brick_mesh = PolygonalMesh::createBrickMesh(Vec3(BL/2.0-0.001, BH/2.0-0.001, 0.1), 10);
	ContactGeometry::TriangleMesh boulder_geo(brick_mesh);

	boulder_body.addContactSurface(Transform(),
        ContactSurface(boulder_geo, material, 0.001));
}


// -----------------Visuals----------------------------------------------------------


void addWallContactsDecor(SimbodyMatterSubsystem& matter)
{
	// Left wall
	const Rotation R_left(Pi, ZAxis);
	matter.Ground().updBody().addDecoration(
		Transform(R_left, Vec3(-BW/2.0, DH/2.0, 0.0)),
        DecorativeBrick(Vec3(BW/2.0, DH/2.0, 0.01)).setColor(Gray));
	// Floor
	const Rotation R_floor(-Pi/2.0, ZAxis);
	matter.Ground().updBody().addDecoration(
		Transform(R_floor, Vec3(DL/2.0, -BW/2.0, 0.0)),
        DecorativeBrick(Vec3(BW/2.0, DL/2.0, 0.01)).setColor(Gray));
	// Ceiling
	const Rotation R_ceiling(Pi/2.0, ZAxis);
	matter.Ground().updBody().addDecoration(
		Transform(R_ceiling, Vec3(DL/2.0, DH + BW/2.0, 0.0)),
        DecorativeBrick(Vec3(BW/2.0, DL/2.0, 0.01)).setColor(Gray));
	// Right wall
	matter.Ground().updBody().addDecoration(Vec3(DL + BW/2.0, DH/2.0, 0.0),
        DecorativeBrick(Vec3(BW/2.0, DH/2.0, 0.01)).setColor(Gray));
}


void addCliffContactDecor(Body& cliff_body) 
{
	cliff_body.addDecoration(Transform(),
        DecorativeBrick(Vec3(cliff_hlx, cliff_hly, 0.01)).setColor(Gray));
}


void addCliffContactDecor(SimbodyMatterSubsystem& matter) 
{
	matter.Ground().updBody().addDecoration(Transform(cliff_center),
        DecorativeBrick(Vec3(cliff_hlx, cliff_hly, 0.01)).setColor(Gray));
}


void addBoulderDecor(Body::Rigid& boulder, const Vec3 color=Green)
{
	boulder.addDecoration(Transform(),
        DecorativeBrick(Vec3(BL/2.0, BH/2.0, 0.01)).setColor(color));
}


class MyReporter : public PeriodicEventReporter {
public:
    MyReporter(const MultibodySystem& system, 
               const CompliantContactSubsystem& complCont,
               Real reportInterval)
    :   PeriodicEventReporter(reportInterval), m_system(system),
        m_compliant(complCont)
    {}

    ~MyReporter() {}

    void handleEvent(const State& state) const override {
        m_system.realize(state, Stage::Position);
        saveEm.push_back(state);
    }
private:
    const MultibodySystem&           m_system;
    const CompliantContactSubsystem& m_compliant;
};