#include "sphinxsys.h"
#include "../../Chrono.h"
#include "../../Log.h"
#include "../../hdf5_IO.h"


using namespace SPH;

// TODO: Add the rest of the walls to chronos
// TODO: Change tank materials
// TODO: Tune contact coefficients

//for geometry
constexpr Real resolution_ref = 5.0e-3;		// < particle spacing
constexpr Real BW = resolution_ref * 4; 	// < boundary width
constexpr Real DL = 0.4;			  		// < tank length
constexpr Real DH = 0.2;				  	// < tank height
constexpr Real DW = 0.1;				  	// < tank width
constexpr Real LL = 0.3 * DL;				// < liquid length
constexpr Real LH = 0.5 * DH;				// < liquid height
constexpr Real LW = DW;				  		// < liquid width
constexpr Real BDL = 2.0e-2;				// < boulder length
constexpr Real BDH = 2.0e-2;				// < boulder height
constexpr Real BDW = 2.0e-2;				// < boulder width
constexpr Real BDx = 0.2 * DL;				// < boulder x position
constexpr Real BDy = 0.5 * DW;				// < boulder y position
constexpr Real BDz = 0.1 * DH;				// < boulder z position

const BoundingBox system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DL + BW, DW + BW, DH + BW));

// Material properties of the fluid
constexpr Real rho0_f = 1000.0;
constexpr Real gravity_g = 9.81;
const Real U_f = 2.0 * sqrt(gravity_g * LH);
const Real c_f = 10.0 * U_f;
constexpr Real mu_f = 0.001;

// Material properties of the solid.
constexpr Real rho0_s = 1000;//2.8e3;								// < Boulder Density [kg/m^3] from paper.
constexpr Real boulder_vol = BDL * BDH * BDW;				// < Boulder Volume [m^3]. (1.5 x 2.0 x 3.0 cm)
constexpr Real boulder_mass = rho0_s * boulder_vol;			// < Boulder Mass [kg].
constexpr Real poisson = 0.3;								// < Poisson's ratio.
constexpr Real Youngs_modulus = 73e9;						// < Young's modulus [Pa].

constexpr int resolution = 0;

constexpr float friction_coef = 0.2f;

std::shared_ptr<ChBody> createChBoulder()
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
	boulder_ch->SetPos(ChVector<>(BDx, BDy, BDz));

	return boulder_ch;
}

std::shared_ptr<ChBody> createChWalls()
{	
	auto tank_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
	tank_mat->SetFriction(friction_coef);

	auto tank_ch = chrono_types::make_shared<ChBody>();
	tank_ch->SetBodyFixed(true);
	tank_ch->SetCollide(true);
	tank_ch->SetDensity(rho0_s);
	tank_ch->GetCollisionModel()->ClearModel();

	// Floor
	tank_ch->GetCollisionModel()->AddBox(tank_mat,
										 0.5 * DL,
										 0.5 * DW,
										 0.5 * BW,
										 ChVector<>(0.5 * DL, 0.5 * DW, -(0.5 * BW)));

	// Top wall
	tank_ch->GetCollisionModel()->AddBox(tank_mat,
										 0.5 * DL,
										 0.5 * DW,
										 0.5 * BW,
										 ChVector<>(0.5 * DL, 0.5 * DW, DH + 0.5*BW));

	// Side wall front
	tank_ch->GetCollisionModel()->AddBox(tank_mat,
										 0.5 * DL,
										 0.5 * BW,
										 0.5 * DH,
										 ChVector<>(0.5 * DL, -(0.5 * BW), 0.5 * DH));
	
	// Side wall back
	tank_ch->GetCollisionModel()->AddBox(tank_mat,
										 0.5 * DL,
										 0.5 * BW,
										 0.5 * DH,
										 ChVector<>(0.5 * DL, DW + 0.5*BW, 0.5 * DH));
	
	// Left wall
	tank_ch->GetCollisionModel()->AddBox(tank_mat,
										 0.5 * BW,
										 0.5 * DW,
										 0.5 * DH,
										 ChVector<>(-(0.5 * BW), 0.5 * DW, 0.5 * DH));
	
	// Right wall
	tank_ch->GetCollisionModel()->AddBox(tank_mat,
										 0.5 * BW,
										 0.5 * DW,
										 0.5 * DH,
										 ChVector<>(DL + 0.5*BW, 0.5 * DW, 0.5 * DH));

	tank_ch->GetCollisionModel()->BuildModel();
	return tank_ch;
}

class WallBoundary : public SolidBodyWithChrono
{
public:
	WallBoundary(SPHSystem &system, const std::string &body_name,
		ChSystemNSC &ch_system)
		: SolidBodyWithChrono(system, body_name, ch_system)
	{
		// tank shape outer shape
		Vec3d outer_dims(DL + 2.0*BW, DW + 2.0*BW, DH + 2.0*BW);
		Vec3d outer_pos(0.5 * DL, 0.5 * DW, 0.5 * DH);
		// inner work
		Vec3d inner_dims(DL, DW, DH);
		Vec3d inner_pos(0.5 * DL, 0.5 * DW, 0.5 * DH);

		body_shape_.add<TriangleMeshShapeBrick>(0.5 * outer_dims, resolution, outer_pos);
		body_shape_.substract<TriangleMeshShapeBrick>(0.5 * inner_dims, resolution, inner_pos);
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
		Vecd wblock_pos(DL - 0.5 * LL, 0.5 * LW, 0.5 * LH);

		body_shape_.add<TriangleMeshShapeBrick>(wblock_dims * 0.5, resolution, wblock_pos);
	}
};

class Boulder : public SolidBodyWithChrono
{
public:
	Boulder(SPHSystem &system, const std::string &body_name,
		ChSystemNSC &ch_system)
		: SolidBodyWithChrono(system, body_name, ch_system)
	{
		Vecd boulder_dims(BDL, BDW, BDH);
		Vecd boulder_pos(BDx, BDy, BDz);

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
