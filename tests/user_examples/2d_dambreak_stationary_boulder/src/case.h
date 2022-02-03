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
Real DL = 2.0; 									/**< Tank length [m]. */
Real DH = 0.4; 									/**< Tank height [m]. */
Real WH = 0.3; 									/**< Water block width [m]. */
Real WL = 0.5; 									/**< Water block height [m]. */
Real BL = 2.0e-2;								/**< Boulder lenght [m]. */
Real BH = 1.5e-2;								/**< Boulder height [m]. */
Real B_x = 0.6;									/**< Boulder initial position x-axis (right edge) [m]. */
Real B_y = 0.04;								/**< Boulder initial position y-axis (right edge) [m]. */
Real particle_spacing_ref = BH / 6.0; 			/**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 6; 			/**< Extending width for BCs. */

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));

/**< Gravity force of fluid. [m/s^2] */
Real gravity_g = 9.81;

/**
 * @brief Material properties of the fluid.
 */
Real rho0_f = 1000.0;								/**< Reference density of fluid [kg/m^3]. */
Real U_f = 2.0 * sqrt(gravity_g * WH);				/**< Characteristic velocity [m/s]. */
Real c_f = 10.0 * U_f;								/**< Reference sound speed [m/s]. */
Real mu_f = 1.0e-3;									/**< Reference dynamic viscocity of fluid [Ns/m^2]. */

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
		std::vector<Vecd> water_block_shape = CreateWaterBlockShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(water_block_shape, ShapeBooleanOps::add);
	}
};


/**
* @brief 	Case dependent material properties definition.
*/
class WaterMaterial : public WeaklyCompressibleFluid
{
public:
	WaterMaterial() : WeaklyCompressibleFluid()
	{
		/** Basic material parameters*/
		rho_0_ = rho0_f;
		c_0_ = c_f;
        mu_ = mu_f;

		/** Compute the derived material parameters*/
		assignDerivedMaterialParameters();
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
		/** Geomtry definition. */
		std::vector<Vecd> outer_shape = CreateOuterWallShape();
		std::vector<Vecd> inner_shape = CreateInnerWallShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(outer_shape, ShapeBooleanOps::add);
		body_shape_->addAPolygon(inner_shape, ShapeBooleanOps::sub);
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
		body_shape_ = new ComplexShape(body_name);
		std::vector<Vecd> boulder_shape = CreateBoulderShape();
		body_shape_->addAPolygon(boulder_shape, ShapeBooleanOps::add);
	}
};
