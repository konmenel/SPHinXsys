/**
 * @brief 	SPHinXsys Library.
 */
#include "sphinxsys.h"

using namespace SPH;

//for geometry
const Real resolution_ref = 5.0e-3;	  	//particle spacing
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
const Real BDx = 0.8;				  	//boulder x position
const Real BDy = DW/2.0;				//boulder y position
const Real BDz = 0.0;				  	//boulder z position

// x="0.9" y="0.24" z="0"
// x="0.12" y="0.12" z="0.45" 


/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DL + BW, DW + BW, DH + BW));

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

/**
 * @brief Contact properties of the simbody.
 */
const Real fK = SimTK::ContactMaterial 							/**< Stiffness Coefficient [Pa] */
		::calcPlaneStrainStiffness(Youngs_modulus, poisson);
const Real fDis = 10.0; // to turn off dissipation
const Real fFac = 0.0;//0.2; // to turn off friction
const Real fVis = 0.0; //0.02; // to turn off viscous friction
const SimTK::ContactMaterial contact_material(fK, fDis, fFac, fFac, fVis);

//	resolution which controls the quality of created polygonalmesh
int resolution(30);

void addSimbodyWallContacts(SimTK::SimbodyMatterSubsystem& matter, 
		const SimTK::ContactCliqueId& clique)
{
	using SimTK::YAxis;
	using SimTK::ZAxis;

	const SimTK::ContactGeometry::HalfSpace half_space;

	// // Left wall
	// const SimTK::Rotation R_left(Pi, ZAxis);
	// matter.Ground().updBody().addContactSurface(
	// 	SimTK::Transform(R_left),
    //     SimTK::ContactSurface(half_space, contact_material)
    //                    .joinClique(clique));
	
	// Floor
	const SimTK::Rotation R_floor(Pi * 0.5, YAxis);
	matter.Ground().updBody().addContactSurface(
		SimTK::Transform(R_floor),
        SimTK::ContactSurface(half_space, contact_material)
                       .joinClique(clique));
	
	// // Ceiling
	// const SimTK::Rotation R_ceiling(Pi * 0.5, YAxis);
	// matter.Ground().updBody().addContactSurface(
	// 	SimTK::Transform(R_ceiling, Vec3d(0.0, 0.0, DH)),
    //     SimTK::ContactSurface(half_space, contact_material)
    //                    .joinClique(clique));
	
	// // Right wall
	// matter.Ground().updBody().addContactSurface(Vec3d(DL, 0.0, 0.0),
    //     SimTK::ContactSurface(half_space, contact_material)
    //                    .joinClique(clique));
	
	// // Front wall
	// const SimTK::Rotation R_front(-Pi * 0.5, ZAxis);
	// matter.Ground().updBody().addContactSurface(
	// 	SimTK::Transform(R_front),
    //     SimTK::ContactSurface(half_space, contact_material)
    //                    .joinClique(clique));
	
	// // Back wall
	// const SimTK::Rotation R_back(Pi * 0.5, ZAxis);
	// matter.Ground().updBody().addContactSurface(
	// 	SimTK::Transform(R_back, Vec3d(0.0, DW, 0.0)),
	// 	SimTK::ContactSurface(half_space, contact_material)
	// 				   .joinClique(clique));
}

void addCliffContactForSimbody(SimTK::SimbodyMatterSubsystem& matter,
		const SimTK::ContactCliqueId& clique) 
{
	Vec3d half_lengths(0.5*(DL - VWx), 0.5 * DW, 0.5 * VWH);
	int resoluton = 0;

	// Create mesh
	SimTK::PolygonalMesh brick_mesh;
	brick_mesh = SimTK::PolygonalMesh::createBrickMesh(half_lengths, resoluton);
	SimTK::ContactGeometry::TriangleMesh cliff_geometry(brick_mesh);

	// Add Contact surface to body
	matter.Ground().updBody().addContactSurface(SimTK::Transform(Vec3d(VWx + 0.5*(DL - VWx), 0.5 * DW, 0.5 * VWH)),
        SimTK::ContactSurface(cliff_geometry, contact_material, resolution_ref)
				.joinClique(clique));
}

void addBoulderContactForSimbody(SimTK::Body::Rigid& boulder_body)
{	
	Vec3d half_lengths(0.5 * BDL, 0.5 * BDW, 0.5 * BDH);
	int resolution = 0;

	// Create mesh
	SimTK::PolygonalMesh brick_mesh;
	brick_mesh = SimTK::PolygonalMesh::createBrickMesh(half_lengths, resolution);
	SimTK::ContactGeometry::TriangleMesh boulder_geo(brick_mesh);

	boulder_body.addContactSurface(SimTK::Transform(),
        SimTK::ContactSurface(boulder_geo, contact_material, resolution_ref));
}

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
		// tank shape
		Vecd halfsize_outer(0.5 * DL + BW, 0.5 * DW + BW, 0.5 * DH + BW);
		Vecd translation_wall(0.5 * DL, 0.5 * DW, 0.5 * DH);
		Vecd halfsize_inner(0.5 * DL, 0.5 * DW, 0.5 * DH);
		body_shape_.add<TriangleMeshShapeBrick>(halfsize_outer, resolution, translation_wall);
		body_shape_.substract<TriangleMeshShapeBrick>(halfsize_inner, resolution, translation_wall);

		// vertical wall
		Vecd halfsize_vwall_outer(0.5*(DL - VWx), 0.5 * DW, 0.5 * VWH);
		Vecd translation_vwall(VWx + 0.5*(DL - VWx), 0.5 * DW, 0.5 * VWH);
		body_shape_.add<TriangleMeshShapeBrick>(halfsize_vwall_outer, resolution, translation_vwall);

		Vecd halfsize_vwall_inner(0.5*(DL - VWx), 0.5*DW + 1.5*BW, 0.5 * VWH);
		translation_vwall = Vecd(VWx + BW + 0.5*(DL - VWx), 0.5 * DW, -BW + 0.5*VWH);
		body_shape_.substract<TriangleMeshShapeBrick>(halfsize_vwall_outer, resolution, translation_vwall);
	}
};

class Boulder : public SolidBody
{
public:
	Boulder(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		Vecd halfsize(0.5 * BDL, 0.5 * BDW, 0.5 * BDH);
		Vecd translation_wall(VWx - BDx + 0.5*BDL, BDy, BDz + 0.5*BDH);
		body_shape_.add<TriangleMeshShapeBrick>(halfsize, resolution, translation_wall);
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
				SimTK::UnitInertia::brick(0.5 * BDL, 0.5 * BDW, 0.5 * BDH)
			);
	}
};