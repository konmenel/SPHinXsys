/**
 * @brief 	SPHinXsys Library.
 */
#include "sphinxsys.h"

using namespace SPH;

//for geometry
const Real resolution_ref = 0.025;	  	//particle spacing
const Real BW = resolution_ref * 4; 	//boundary width
const Real DL = 1.6;			  		//tank length
const Real DH = 0.6;				  	//tank height
const Real DW = 0.8;				  	//tank width
const Real LL = 0.4;				  	//liquid length
const Real LH = 0.4;				  	//liquid height
const Real LW = 0.8;				  	//liquid width

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DL + BW, DW + BW, DH + BW));

//for material properties of the fluid
const Real rho0_f = 1000.0;
const Real gravity_g = 9.81;
const Real U_f = 2.0 * sqrt(gravity_g * LH);
const Real c_f = 10.0 * U_f;
const Real mu_f = 0.001;

//	resolution which controls the quality of created polygonalmesh
int resolution(50);

//	define the fluid body
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &system, const std::string &body_name)
		: FluidBody(system, body_name)
	{
		Vecd halfsize_water(0.5 * LL, 0.5 * LW, 0.5 * LH);
		Vecd translation_water = halfsize_water;

		body_shape_.add<TriangleMeshShapeBrick>(halfsize_water, resolution, translation_water);
	}
};
//	define the static solid wall boundary
class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		Vecd halfsize_outer(0.5 * DL + BW, 0.5 * DW + BW, 0.5 * DH + BW);
		Vecd translation_wall(0.5 * DL, 0.5 * DW, 0.5 * DH);
		Vecd halfsize_inner(0.5 * DL, 0.5 * DW, 0.5 * DH);
		body_shape_.add<TriangleMeshShapeBrick>(halfsize_outer, resolution, translation_wall);
		body_shape_.substract<TriangleMeshShapeBrick>(halfsize_inner, resolution, translation_wall);
	}
};
