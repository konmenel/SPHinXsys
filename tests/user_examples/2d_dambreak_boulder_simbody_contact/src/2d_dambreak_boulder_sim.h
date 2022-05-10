 /**
  * @brief 	SPHinXsys Library.
  */
#pragma once
#include "sphinxsys.h"
  /**
 * @brief Namespace cite here.
 */
using namespace SPH;

// TODO: add the vertical wall
// TODO: add dumping area at the end
// TODO: add observers, wave probes etc.

/**
 * @brief Basic geometry parameters and numerical setup.
 */
const Real DL = 3.0; 								/**< Tank length [m]. */
const Real DH = 0.5; 								/**< Tank height [m]. */
const Real WH = 0.25; 								/**< Water block width [m]. */
const Real WL = 0.5; 								/**< Water block height [m]. */
const Real wall_H = 0.2;							/**< Verical wall height [m] */
const Real wall_position = 2.0;						/**< Position of the verical wall [m] (x direction) */ 
const Real BL = 2.0e-2;								/**< Boulder lenght [m]. */
const Real BH = 1.5e-2;								/**< Boulder height [m]. */
const Real B_x = wall_position - 0.1;				/**< Boulder initial position x-axis (right edge) [m]. */
const Real B_y = 0.1;								/**< Boulder initial position y-axis (bottom edge) [m]. */
const Real particle_spacing_ref = BH / 5.0; 		/**< Initial reference particle spacing. */
const Real BW = particle_spacing_ref * 4.0; 		/**< Extending width for BCs. */

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));

const Real gravity_g = 9.81;		/**< Gravity force of fluid. [m/s^2] */

/**
 * @brief Material properties of the fluid.
 */
const Real rho0_f = 1000.0;								/**< Reference density of fluid [kg/m^3]. */
const Real U_f = 2.0 * sqrt(gravity_g * WH);			/**< Characteristic velocity [m/s]. */
const Real c_f = 10.0 * U_f;							/**< Reference sound speed [m/s]. */
const Real mu_f = 1e-3;									/**< Reference dynamic viscocity of fluid [Ns/m^2]. */

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
const Vec3d cliff_center(wall_position + cliff_hlx, cliff_hly, 0.0);

/**
 * @brief Contact properties of the simbody.
 */
const Real fK = SimTK::ContactMaterial 							/**< Stiffness Coefficient [Pa] */
		::calcPlaneStrainStiffness(Youngs_modulus, poisson);
const Real fDis = 10.0; // to turn off dissipation
const Real fFac = 0.2; // to turn off friction
const Real fVis = 0.02; // to turn off viscous friction
const SimTK::ContactMaterial contact_material(fK, fDis, fFac, fFac, fVis);
const Real surface_thickness = 1.0;
const Real hlz = 0.01;					/**< Brick half-length in z direction (Needed for Simbody) */					


void addSimbodyWallContacts(SimTK::SimbodyMatterSubsystem& matter, 
		const SimTK::ContactCliqueId& clique)
{
	using SimTK::ZAxis;

	const SimTK::ContactGeometry::HalfSpace half_space;
	// Left wall
	const SimTK::Rotation R_left(Pi, ZAxis);
	matter.Ground().updBody().addContactSurface(
		SimTK::Transform(R_left),
        SimTK::ContactSurface(half_space, contact_material)
                       .joinClique(clique));
	// Floor
	const SimTK::Rotation R_floor(-Pi/2.0, ZAxis);
	matter.Ground().updBody().addContactSurface(
		SimTK::Transform(R_floor),
        SimTK::ContactSurface(half_space, contact_material)
                       .joinClique(clique));
	// Ceiling
	const SimTK::Rotation R_ceiling(Pi/2.0, ZAxis);
	matter.Ground().updBody().addContactSurface(
		SimTK::Transform(R_ceiling, Vec3d(0.0, DH, 0.0)),
        SimTK::ContactSurface(half_space, contact_material)
                       .joinClique(clique));
	// Right wall
	matter.Ground().updBody().addContactSurface(Vec3d(DL, 0.0, 0.0),
        SimTK::ContactSurface(half_space, contact_material)
                       .joinClique(clique));
}


void addCliffContactForSimbody(SimTK::SimbodyMatterSubsystem& matter,
		const SimTK::ContactCliqueId& clique) 
{
	Vec3d half_lengths(cliff_hlx, cliff_hly, hlz);
	int resoluton = 0;

	// Create mesh
	SimTK::PolygonalMesh brick_mesh;
	brick_mesh = SimTK::PolygonalMesh::createBrickMesh(half_lengths, resoluton);
	SimTK::ContactGeometry::TriangleMesh cliff_geometry(brick_mesh);

	// Add Contact surface to body
	matter.Ground().updBody().addContactSurface(SimTK::Transform(cliff_center),
        SimTK::ContactSurface(cliff_geometry, contact_material, surface_thickness)
				.joinClique(clique));
}


void addBoulderContactForSimbody(SimTK::Body::Rigid& boulder_body)
{	
	Vec3d half_lengths(BL/2.0, BH/2.0, hlz);
	int resolution = 0;

	// Create mesh
	SimTK::PolygonalMesh brick_mesh;
	brick_mesh = SimTK::PolygonalMesh::createBrickMesh(half_lengths, resolution);
	SimTK::ContactGeometry::TriangleMesh boulder_geo(brick_mesh);

	boulder_body.addContactSurface(SimTK::Transform(),
        SimTK::ContactSurface(boulder_geo, contact_material, surface_thickness));
}


/**
* @brief 	Create a water block shape.
*/
std::vector<Vecd> createWaterBlockShape()
{
	//geometry
	std::vector<Vecd> water_block_shape;
	water_block_shape.push_back(Vecd(0.0, 0.0));
	water_block_shape.push_back(Vecd(0.0, WH));
	water_block_shape.push_back(Vecd(WL, WH));
	water_block_shape.push_back(Vecd(WL, 0.0));
	water_block_shape.push_back(Vecd(0.0, 0.0));
	return water_block_shape;
}

/** 
* @brief	Create outer wall shape.
*/
std::vector<Vecd> createOuterWallShape()
{
	std::vector<Vecd> outer_wall_shape;
	outer_wall_shape.push_back(Vecd(-BW, -BW));
	outer_wall_shape.push_back(Vecd(-BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
	outer_wall_shape.push_back(Vecd(DL + BW, -BW));
	outer_wall_shape.push_back(Vecd(-BW, -BW));
	return outer_wall_shape;
}

/**
* @brief 	Create inner wall shape.
*/
std::vector<Vecd> createInnerWallShape()
{
	std::vector<Vecd> inner_wall_shape;
	inner_wall_shape.push_back(Vecd(0.0, 0.0));
	inner_wall_shape.push_back(Vecd(0.0, DH));
	inner_wall_shape.push_back(Vecd(DL, DH));
	inner_wall_shape.push_back(Vecd(DL, 0.0));
	inner_wall_shape.push_back(Vecd(0.0, 0.0));
	return inner_wall_shape;
}
/**
* @brief 	Create verical wall shape.
*/
std::vector<Vecd> createVerWallShape()
{
	std::vector<Vecd> inner_wall_shape;
	inner_wall_shape.push_back(Vecd(wall_position, 0.0));
	inner_wall_shape.push_back(Vecd(wall_position, wall_H));
	inner_wall_shape.push_back(Vecd(DL, wall_H));
	inner_wall_shape.push_back(Vecd(DL, 0.0));
	inner_wall_shape.push_back(Vecd(wall_position, 0.0));
	return inner_wall_shape;
}

/**
* @brief 	Create boulder shape.
*/
std::vector<Vecd> createBoulderShape()
{
	std::vector<Vecd> boulder;
	boulder.push_back(Vecd(B_x - BL, B_y));
	boulder.push_back(Vecd(B_x - BL, B_y + BH));
	boulder.push_back(Vecd(B_x, B_y + BH));
	boulder.push_back(Vecd(B_x, B_y));
	boulder.push_back(Vecd(B_x - BL, B_y));
	return boulder;
}

MultiPolygon createBoulderSimbodyConstrainShape()
{
	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(createBoulderShape(), ShapeBooleanOps::add);
	return multi_polygon;
};

/**
* @brief 	Fluid body definition.
*/
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem& sph_system, string body_name)
		: FluidBody(sph_system, body_name)
	{
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};

/**
 * @brief 	Wall boundary body definition.
 */
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &sph_system, string body_name)
		: SolidBody(sph_system, body_name)
	{
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createOuterWallShape(), ShapeBooleanOps::add);
		multi_polygon.addAPolygon(createInnerWallShape(), ShapeBooleanOps::sub);
		multi_polygon.addAPolygon(createVerWallShape(), ShapeBooleanOps::add);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};


/**
 * @brief 	Boulder body definition.
 */
class Boulder : public SolidBody
{
public:
	Boulder(SPHSystem &sph_system, string body_name)
		: SolidBody(sph_system, body_name)
	{
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createBoulderShape(), ShapeBooleanOps::add);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};

/**
* @brief 	Create boulder body for simbody
*/
class BoulderSystemForSimbody : public SolidBodyPartForSimbody
{
public:
	BoulderSystemForSimbody(SolidBody &solid_body,
						 	const std::string &constrained_region_name,
							Shape& shape)
		: SolidBodyPartForSimbody(solid_body, constrained_region_name, shape)
	{
		body_part_mass_properties_ = mass_properties_ptr_keeper_
			.createPtr<SimTK::MassProperties>(
				boulder_mass, 
				SimTK::Vec3(0.0), 
				SimTK::UnitInertia::brick(0.5 * BL, 0.5 * BH, 0.0)
			);
	}
};
