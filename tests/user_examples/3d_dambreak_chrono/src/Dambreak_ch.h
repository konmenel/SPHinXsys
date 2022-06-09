#include "sphinxsys.h"
#include "Chrono.h"
#include "Log.h"

using namespace SPH;

// TODO: Add the rest of the walls to chronos
// TODO: True to create each wall of the tank seperately
// TODO: Change tank materials
// TODO: Tune contact coefficients

//for geometry
const Real resolution_ref = 3.75e-3;	//particle spacing
const Real BW = resolution_ref * 4; 	//boundary width
const Real DL = 4.0;			  		//tank length
const Real DH = 0.5;				  	//tank height
const Real DW = 0.3;				  	//tank width
const Real LL = 1.5;				  	//liquid length
const Real LH = 0.2;				  	//liquid height
const Real LW = 0.3;				  	//liquid width
const Real BDL = 3.0e-2;				//boulder length
const Real BDH = 1.5e-2;				//boulder height
const Real BDW = 2.0e-2;				//boulder width
const Real VWx = 2.0;					//verical wall x position
const Real VWH = 0.2;					//verical wall height
const Real BDx = 0.08;				  	//boulder x position
const Real BDy = DW * 0.5;				//boulder y position
const Real BDz = 0.05;				  	//boulder z position

const BoundingBox system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DL + BW, DW + BW, DH + BW));

// Material properties of the fluid
const Real rho0_f = 1000.0;
const Real gravity_g = 9.81;
const Real U_f = 2.0 * sqrt(gravity_g * LH);
const Real c_f = 10.0 * U_f;
const Real mu_f = 0.001;

// Material properties of the solid.
const Real rho0_s = 2.8e3;								/**< Boulder Density [kg/m^3] from paper. */
const Real boulder_vol = BDL * BDH * BDW;				/**< Boulder Volume [m^3]. (1.5 x 2.0 x 3.0 cm) */
const Real boulder_mass = rho0_s * boulder_vol;			/**< Boulder Mass [kg]. */
const Real poisson = 0.3;								/**< Poisson's ratio. */
const Real Youngs_modulus = 73e9;						/**< Young's modulus [Pa]. */

const int resolution = 0;

const float friction_coef = 0.2f;
const auto collision_type = chrono::collision::ChCollisionSystemType::BULLET;

std::shared_ptr<ChBody> addBoulderCh(ChSystem &ch_system)
{
	auto boulder_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
	boulder_mat->SetFriction(friction_coef);

	auto boulder_ch = chrono_types::make_shared<ChBodyEasyBox>(	BDL,
															BDW,
															BDH,
															rho0_s,
															false,
															true,
															boulder_mat);
	boulder_ch->SetPos(ChVector<>(VWx - BDx - 0.5*BDL, BDy, BDz + 0.5*BDH));

	ch_system.AddBody(boulder_ch);
	return boulder_ch;
}

void addWallsCh(ChSystem &ch_system)
{	
	auto tank_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
	tank_mat->SetFriction(friction_coef);

	auto floor1_ch = chrono_types::make_shared<ChBodyEasyBox>(	VWx,
																DW,
																BW,
																rho0_s,
																false,
																true,
																tank_mat);
	floor1_ch->SetPos(ChVector<>(0.5 * VWx, 0.5 * DW, -(0.5 * BW)));
	floor1_ch->SetBodyFixed(true);
	ch_system.AddBody(floor1_ch);


	auto floor2_ch = chrono_types::make_shared<ChBodyEasyBox>(	DL - VWx,
																DW,
																BW,
																rho0_s,
																false,
																true,
																tank_mat);
	floor2_ch->SetPos(ChVector<>(0.5*(DL + VWx), 0.5 * DW, VWH - 0.5*BW));
	floor2_ch->SetBodyFixed(true);
	ch_system.AddBody(floor2_ch);

	auto verical_wall_ch = chrono_types
		::make_shared<ChBodyEasyBox>(	BW,
										DW,
										VWH - BW,
										rho0_s,
										false,
										true,
										tank_mat);
	verical_wall_ch->SetPos(ChVector<>(VWx + 0.5*BW, 0.5 * DW, 0.5*(VWH - BW)));
	verical_wall_ch->SetBodyFixed(true);
	ch_system.AddBody(verical_wall_ch);
	return;
}

class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		// // tank shape outer shape
		Vec3d outer_dims(DL + 2.0*BW, DW + 2.0*BW, DH + 2.0*BW);
		Vec3d outer_pos(0.5 * DL, 0.5 * DW, 0.5 * DH);
		// // Remove parts at the right part of the cliff
		Vec3d linner_dims(VWx, DW, DH);
		Vec3d linner_pos(0.5 * VWx, 0.5 * DW, 0.5 * DH);
		// // Remove parts above cliff
		Vec3d rinner_dims(DL - VWx, DW, DH - VWH);
		Vec3d rinner_pos(0.5*(DL + VWx), 0.5*DW, VWH + 0.5*(DH - VWH));
		// remove inside of cliff
		Vec3d wallinner_dims(DL - VWx, DW + 2.0*BW, VWH);
		Vec3d wallinner_pos(0.5*(DL + VWx) + BW , 0.5 * DW, 0.5*VWH - BW);

		body_shape_.add<TriangleMeshShapeBrick>(0.5 * outer_dims, resolution, outer_pos);
		body_shape_.substract<TriangleMeshShapeBrick>(0.5 * linner_dims, resolution, linner_pos);
		body_shape_.substract<TriangleMeshShapeBrick>(0.5 * rinner_dims, resolution, rinner_pos);
		body_shape_.substract<TriangleMeshShapeBrick>(0.5 * wallinner_dims, resolution, wallinner_pos);
	}
};

// Define the fluid body
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &system, const std::string &body_name)
		: FluidBody(system, body_name)
	{
		Vecd wblock_dims(LL, LW, LH);
		Vecd wblock_pos(0.5 * LL, 0.5 * LW, 0.5 * LH);

		body_shape_.add<TriangleMeshShapeBrick>(wblock_dims * 0.5, resolution, wblock_pos);
	}
};

class Boulder : public SolidBody
{
public:
	Boulder(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		Vecd boulder_dims(BDL, BDW, BDH);
		Vecd boulder_pos(VWx - BDx - 0.5*BDL, BDy, BDz + 0.5*BDH);

		body_shape_.add<TriangleMeshShapeBrick>(boulder_dims * 0.5, resolution, boulder_pos);
	}
};

/**
* @brief 	Create boulder body for Chrono
*/
class BoulderSystemForChrono : public BodyRegionByParticle
{
public:
	BoulderSystemForChrono(	SolidBody &solid_body,
						 	const std::string &constrained_region_name,
							Shape& shape)
		: BodyRegionByParticle(solid_body, constrained_region_name, shape) {}

	virtual ~BoulderSystemForChrono() {}
};
