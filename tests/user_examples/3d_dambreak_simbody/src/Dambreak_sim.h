/**
 * @brief 	SPHinXsys Library.
 */
#include "sphinxsys.h"

using namespace SPH;

//for geometry
const Real resolution_ref = 3.75e-3;	  	//particle spacing
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
const Real surface_thickness = 1.0;

/**
 * @brief Contact properties of the simbody.
 */
const Real fK = SimTK::ContactMaterial 							/**< Stiffness Coefficient [Pa] */
		::calcPlaneStrainStiffness(Youngs_modulus, poisson);
const Real fDis = 10.0; // to turn off dissipation
const Real fFac = 0.3; // to turn off friction
const Real fVis = 0.0; //0.02; // to turn off viscous friction
const SimTK::ContactMaterial contact_material(fK, fDis, fFac, fFac, fVis);

//	resolution which controls the quality of created polygonalmesh
const int resolution = 0;

static void createBrickMesh(const Vecd hf, SimTK::PolygonalMesh &brick)
{
	SimTK::Array_<Vecd> vertices;
	vertices.push_back(Vecd( hf[0],  hf[1],  hf[2]));
	vertices.push_back(Vecd( hf[0],  hf[1], -hf[2]));
	vertices.push_back(Vecd( hf[0], -hf[1],  hf[2]));
	vertices.push_back(Vecd( hf[0], -hf[1], -hf[2]));
	vertices.push_back(Vecd(-hf[0],  hf[1],  hf[2]));
	vertices.push_back(Vecd(-hf[0],  hf[1], -hf[2]));
	vertices.push_back(Vecd(-hf[0], -hf[1],  hf[2]));
	vertices.push_back(Vecd(-hf[0], -hf[1], -hf[2]));

	SimTK::Array_<int> facesIndeces = {
		0,2,3,1,
		1,5,4,0,
		0,4,6,2,
		2,6,7,3,
		3,7,5,1,
		4,5,7,6
	};
	
	for (size_t i = 0; i < vertices.size(); ++i) {
		brick.addVertex(vertices[i]);
	}
	for (size_t i = 0; i < facesIndeces.size(); i += 4) {
		const SimTK::Array_<int> verts(&facesIndeces[i], &facesIndeces[i]+4);
		brick.addFace(verts);
	}
}

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
	// int resolution = 0;

	// Create mesh
	SimTK::PolygonalMesh brick_mesh;
	// brick_mesh = SimTK::PolygonalMesh::createBrickMesh(half_lengths, resolution);
	createBrickMesh(half_lengths, brick_mesh);
	SimTK::ContactGeometry::TriangleMesh cliff_geometry(brick_mesh);

	// Add Contact surface to body
	matter.Ground().updBody().addContactSurface(SimTK::Transform(Vec3d(VWx + 0.5*(DL - VWx), 0.5 * DW, 0.5 * VWH)),
        SimTK::ContactSurface(cliff_geometry, contact_material, surface_thickness)
				.joinClique(clique));
}

void addBoulderContactForSimbody(SimTK::Body::Rigid& boulder_body)
{	
	Vec3d half_lengths(0.5 * BDL, 0.5 * BDW, 0.5 * BDH);
	// int resolution = 1;

	// Create mesh
	SimTK::PolygonalMesh brick_mesh;
	// brick_mesh = SimTK::PolygonalMesh::createBrickMesh(half_lengths, resolution);
	createBrickMesh(half_lengths, brick_mesh);
	SimTK::ContactGeometry::TriangleMesh boulder_geo(brick_mesh);
	// SimTK::ContactGeometry::Brick boulder_geo(half_lengths);

	boulder_body.addContactSurface(SimTK::Transform(),
        SimTK::ContactSurface(boulder_geo, contact_material, surface_thickness));
}

//	define the fluid body
class WaterBlock : public FluidBody
{
public:
	WaterBlock(SPHSystem &system, const std::string &body_name)
		: FluidBody(system, body_name)
	{
		Vecd halfsize_water(0.5 * LL, 0.5 * LW, 0.5 * LH);
		Vecd translation_water(0.5 * LL, 0.5 * LW, 0.5 * LH);

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

class Boulder : public SolidBody
{
public:
	Boulder(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		Vecd halfsize(0.5 * BDL, 0.5 * BDW, 0.5 * BDH);
		Vecd translation_wall(VWx - BDx - 0.5*BDL, BDy, BDz + 0.5*BDH);
		
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

/**
 * @brief Class to logs output to a file and stdout at the same time
 * 
 */
class LogOutput
{
public:
	LogOutput(const std::string &file_name)
	{
		log_file_.open(file_name);
		is_opened_ = true;
	}

	~LogOutput()
	{
		if (is_opened_)
		{
			log_file_.close();
		}
	}

	void close()
	{
		if (is_opened_)
		{
			log_file_.close();
			is_opened_ = false;
		}
	}

	// Overload the << operator to write to both stdout and log file
	template<typename T>
	LogOutput& operator<<(const T& t)
	{
		std::cout << t;
		log_file_ << t;

		return *this;
	}

	// Overload the << operator to handle std::endl
	LogOutput& operator<<(std::ostream& (*f)(std::ostream&))
	{
		f(std::cout);
		f(log_file_);

		return *this;
	}
	
private:
	bool is_opened_;
	std::ofstream log_file_;
};
