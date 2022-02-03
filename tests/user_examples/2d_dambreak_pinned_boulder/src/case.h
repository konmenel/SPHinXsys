 /**
  * @brief 	SPHinXsys Library.
  */
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
Real DL = 3.0; 									/**< Tank length [m]. */
Real DH = 0.4; 									/**< Tank height [m]. */
Real WH = 0.2; 									/**< Water block width [m]. */
Real WL = 0.5; 									/**< Water block height [m]. */
Real BL = 1.5e-2;								/**< Boulder lenght [m]. */
Real BH = 2.0e-2;								/**< Boulder height [m]. */
Real B_x = 0.6;									/**< Boulder initial position x-axis (right edge) [m]. */
Real B_y = 0.05;								/**< Boulder initial position y-axis (bottom edge) [m]. */
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
 * @brief Material properties of the solid.
 */
Real rho0_s = 2.8e3;								/**< Boulder Density [kg/m^3]. */
Real boulder_vol  = 9e-6;							/**< Boulder Volume [m^3]. (1.5 x 2.0 x 3.0 cm) */
Real boulder_mass = rho0_s * boulder_vol;			/**< Boulder Mass [kg]. */
Real poisson = 0.25;								/**< Poisson's ratio. */
Real Youngs_modulus = 73e9;						//**< Young's modulus [Pa]. */

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
* @brief create a damping zone
*/
// std::vector<Vecd> CreatDampingBufferShape()
// {
// 	std::vector<Vecd> pnts;
// 	pnts.push_back(Vecd(DL - 5.0,  0.356 - BW));
//     pnts.push_back(Vecd(DL - 5.0,  DH));
//     pnts.push_back(Vecd(DL + BW,   DH));
//     pnts.push_back(Vecd(DL + BW,   0.356 - BW));
//     pnts.push_back(Vecd(DL - 5.0,  0.356 - BW));

// 	return pnts;
// }

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

/**
 * @brief Define boulder material.
 */
class BoulderMaterial : public LinearElasticSolid
{
public:
	BoulderMaterial() : LinearElasticSolid()
	{
		rho_0_ = rho0_s;
		E_0_ = Youngs_modulus;
		nu_ = poisson;

		assignDerivedMaterialParameters();
	}
};

/**
* @brief 	Create fish head for constraint
*/
class BoulderSystemForSimbody : public SolidBodyPartForSimbody
{
	void tagBodyPart() override
	{
		BodyPartByParticle::tagBodyPart();
		//Vecd mass_center = Vecd(7.92, 0.355); // 0.3355
		//initial_mass_center_ = SimTK::Vec3(mass_center[0], mass_center[1], 0.0);
		/** UnitInertia_ (const RealP &xx, const RealP &yy, const RealP &zz)
		 * 	Create a principal unit inertia matrix (only non-zero on diagonal). 
		 */
		// Real Iz = (BH * BH + BL * BL) / 12.0; 	/**< Unit Inertia of boulder divided by its mass. */
		body_part_mass_properties_
			= new SimTK::MassProperties(
				boulder_mass, 
				SimTK::Vec3(0.0), 
				SimTK::UnitInertia::brick(BL/2.0, BH/2.0, 0.0)
				);
	}
public:
	BoulderSystemForSimbody(SolidBody* solid_body,
		string constrained_region_name, Real solid_body_density)
		: SolidBodyPartForSimbody(solid_body, constrained_region_name)
	{
		body_part_shape_ = new ComplexShape(constrained_region_name);
		std::vector<Vecd> boulder_shape = CreateBoulderShape();
		body_part_shape_->addAPolygon(boulder_shape, ShapeBooleanOps::add);
		/** tag the constrained particle. */
		tagBodyPart();
	}
};

/**
 * @brief 	Fluid observer body definition.
 */
// class FluidObserver : public FictitiousBody
// {
// public:
// 	FluidObserver(SPHSystem &sph_system, string body_name)
// 		: FictitiousBody(sph_system, body_name)
// 	{
// 		body_input_points_volumes_.push_back(make_pair(Vecd(DL, 0.2), 0.0));
// 	}
// };