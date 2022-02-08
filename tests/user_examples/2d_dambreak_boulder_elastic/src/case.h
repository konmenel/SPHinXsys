 /**
  * @brief 	SPHinXsys Library.
  */
#include "sphinxsys.h"
  /**
 * @brief Namespace cite here.
 */
using namespace SPH;


/**
 * @brief Basic geometry parameters and numerical setup.
 */
const Real DL = 3.5; 								/**< Tank length [m]. */
const Real DH = 0.5; 								/**< Tank height [m]. */
const Real WH = 0.2; 								/**< Water block width [m]. */
const Real WL = 0.5; 								/**< Water block height [m]. */
const Real wall_H = 0.2;							/**< Verical wall height [m] */
const Real wall_position = 2.0;						/**< Position of the verical wall [m] (x direction) */ 
const Real BL = 2.0e-2;								/**< Boulder lenght [m]. */
const Real BH = 1.5e-2;								/**< Boulder height [m]. */
const Real B_x = wall_position - 0.1;				/**< Boulder initial position x-axis (right edge) [m]. */
const Real B_y = 0.01;								/**< Boulder initial position y-axis (right edge) [m]. */
const Real particle_spacing_ref = BH / 6.0; 		/**< Initial reference particle spacing. */
const Real BW = particle_spacing_ref * 6.0; 		/**< Extending width for BCs. */

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));

/**< Gravity force of fluid. [m/s^2] */
const Real gravity_g = 9.81;

/**
 * @brief Material properties of the fluid.
 */
const Real rho0_f = 1000.0;								/**< Reference density of fluid [kg/m^3]. */
const Real U_f = 2.0 * sqrt(gravity_g * WH);			/**< Characteristic velocity [m/s]. */
const Real c_f = 10.0 * U_f;							/**< Reference sound speed [m/s]. */
const Real mu_f = 1.0e-3;								/**< Reference dynamic viscocity of fluid [Ns/m^2]. */

/**
 * @brief Material properties of the solid.
 */
const Real rho0_s = 2.8e3;								/**< Boulder Density [kg/m^3] from paper. */
const Real boulder_vol = BL * BH;						/**< Boulder Volume [m^3]. (1.5 x 2.0 x 3.0 cm) */
const Real boulder_mass = rho0_s * boulder_vol;			/**< Boulder Mass [kg]. */
const Real poisson = 0.25;								/**< Poisson's ratio. */
const Real Youngs_modulus = 73e5;						/**< Young's modulus [Pa]. */
const Real physical_viscosity = 100000.0;

/**
* @brief 	Create a water block shape.
*/
std::vector<Vecd> CreateWaterBlockShape()
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
* @brief 	Create verical wall shape.
*/
std::vector<Vecd> CreateVerWallShape()
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
* @brief 	Fluid body definition.
*/
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem& sph_system, string body_name)
		: FluidBody(sph_system, body_name)
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		std::vector<Vecd> water_block_shape = CreateWaterBlockShape();
		multi_polygon.addAPolygon(water_block_shape, ShapeBooleanOps::add);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};


/**
* @brief 	Case dependent material properties definition.
*/
class WaterMaterial : public WeaklyCompressibleFluid
{
public:
		WaterMaterial() : WeaklyCompressibleFluid(rho0_f, c_f, mu_) {}
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
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		std::vector<Vecd> outer_shape = CreateOuterWallShape();
		std::vector<Vecd> inner_shape = CreateInnerWallShape();
		std::vector<Vecd> verical_wall_shape = CreateVerWallShape();
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
