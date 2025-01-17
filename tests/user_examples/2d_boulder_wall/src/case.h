 /**
  * @brief 	SPHinXsys Library.
  */
#include "sphinxsys.h"
  /**
 * @brief Namespace cite here.
 */
using namespace SPH;
using std::cout;
using std::endl;


/**
 * @brief Basic geometry parameters and numerical setup.
 */
const Real DL = 0.4; 								/**< Tank length [m]. */
const Real DH = 0.5; 								/**< Tank height [m]. */
const Real BL = 2.0e-2;								/**< Boulder lenght [m]. */
const Real BH = 1.5e-2;								/**< Boulder height [m]. */
const Real B_x = DL / 2.0;							/**< Boulder initial position x-axis (right edge) [m]. */
const Real B_y = 0.1;								/**< Boulder initial position y-axis (right edge) [m]. */
const Real particle_spacing_ref = BH / 24.0; 		/**< Initial reference particle spacing. */
const Real BW = particle_spacing_ref * 12.0; 		/**< Extending width for BCs. */

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));

/**< Gravity force of fluid. [m/s^2] */
const Real gravity_g = 9.81;

/**
 * @brief Material properties of the solid.
 */
const Real rho0_s = 2.8e3;								/**< Boulder Density [kg/m^3] from paper. */
const Real boulder_vol = BL * BH;						/**< Boulder Volume [m^3]. (1.5 x 2.0 x 3.0 cm) */
const Real boulder_mass = rho0_s * boulder_vol;			/**< Boulder Mass [kg]. */
const Real poisson = 0.4;								/**< Poisson's ratio. */
const Real Youngs_modulus = 1e9;						/**< Young's modulus [Pa]. */
const Real physical_viscosity = 1e5;

/** 
* @brief	Create outer wall shape.
*/
std::vector<Vecd> CreateOuterWallShape()
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
std::vector<Vecd> CreateInnerWallShape()
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
* @brief 	Create boulder shape.
*/
std::vector<Vecd> CreateBoulderShape()
{
	std::vector<Vecd> boulder;
	boulder.push_back(Vecd(B_x - BL, B_y));
	boulder.push_back(Vecd(B_x - BL, B_y + BH));
	boulder.push_back(Vecd(B_x, B_y + BH));
	boulder.push_back(Vecd(B_x, B_y));
	boulder.push_back(Vecd(B_x - BL, B_y));
	return boulder;
}

/**
 * @brief 	Wall boundary body definition.
 */
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &sph_system, string body_name)
		: SolidBody(sph_system, body_name)
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		std::vector<Vecd> outer_shape = CreateOuterWallShape();
		std::vector<Vecd> inner_shape = CreateInnerWallShape();
		multi_polygon.addAPolygon(outer_shape, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(inner_shape, ShapeBooleanOps::sub);
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
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		std::vector<Vecd> boulder_shape = CreateBoulderShape();
		multi_polygon.addAPolygon(boulder_shape, ShapeBooleanOps::add);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};

/**
 * @brief Define boulder material.
 */
class BoulderMaterial : public LinearElasticSolid
{
public:
	BoulderMaterial() : LinearElasticSolid(rho0_s, Youngs_modulus, poisson) {}
};

MultiPolygon CreateBoulderMultiShape()
{
	MultiPolygon multi_polygon;
	std::vector<Vecd> boulder_shape = CreateBoulderShape();
	multi_polygon.addAPolygon(boulder_shape, ShapeBooleanOps::add);
	return multi_polygon;
}

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
				SimTK::UnitInertia::brick(BL/2.0, BH/2.0, 0.0)
			);
	}
};
