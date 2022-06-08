#include "sphinxsys.h"
#include "Chrono.h"

using namespace SPH;

// TODO: Add the rest of the walls to chronos
// TODO: True to create each wall of the tank seperately
// TODO: Change tank materials
// TODO: Tune contact coefficients

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

const int resolution = 20;

auto collition_model = chrono_types::make_shared<chrono::collision::ChCollisionModelBullet>();
const float friction_coef = 0.2f;

std::shared_ptr<ChBody> addBoulderCh(ChSystem &ch_system)
{
	auto box_mat = chrono_types::make_shared<SurfMaterialCh>();
	box_mat->SetFriction(friction_coef);
#if SMOOTH_CONTACT
	box_mat->SetYoungModulus(Youngs_modulus);
	box_mat->SetPoissonRatio(poisson);
#endif // end if SMOOTH_CONTACT

	auto box_ch = chrono_types::make_shared<ChBodyEasyBox>(	BDL,
															BDW,
															BDH,
															rho0_s,
															false,
															true,
															box_mat,
															collition_model);
	box_ch->SetName("Box");
	box_ch->SetPos(ChVector<>(BDx, BDy, BDz));

	ch_system.AddBody(box_ch);
	return box_ch;
}

void addWallsCh(ChSystem &ch_system)
{	
	auto tank_mat = chrono_types::make_shared<SurfMaterialCh>();
	tank_mat->SetFriction(friction_coef);
#if SMOOTH_CONTACT
	tank_mat->SetYoungModulus(Youngs_modulus);
	tank_mat->SetPoissonRatio(poisson);
#endif // end if SMOOTH_CONTACT

	auto floor1_ch = chrono_types::make_shared<ChBodyEasyBox>(	VWx,
																DW,
																BW,
																rho0_s,
																false,
																true,
																tank_mat,
																collition_model);
	floor1_ch->SetName("Floor1");
	floor1_ch->SetPos(ChVector<>(VWx * 0.5, DW * 0.5, -BW * 0.5));
	floor1_ch->SetBodyFixed(true);
	ch_system.AddBody(floor1_ch);


	auto floor2_ch = chrono_types::make_shared<ChBodyEasyBox>(	DL - VWx,
																DW,
																BW,
																rho0_s,
																false,
																true,
																tank_mat,
																collition_model);
	floor2_ch->SetName("Floor2");
	floor2_ch->SetPos(ChVector<>((DL - VWx)*0.5, DW * 0.5, VWH - BW*0.5));
	floor2_ch->SetBodyFixed(true);
	ch_system.AddBody(floor2_ch);

	auto verical_wall_ch = chrono_types
		::make_shared<ChBodyEasyBox>(	BW,
										DW,
										VWH,
										rho0_s,
										false,
										true,
										tank_mat,
										collition_model);
	verical_wall_ch->SetName("Cliff");
	verical_wall_ch->SetPos(ChVector<>(VWx + BW*0.5, DW * 0.5, VWH*0.5 - BW));
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
		// tank shape outer shape
		Vecd halfsize_outer(0.5*DL + BW, 0.5*DW + BW, 0.5*DH + BW);
		Vecd translation_outer(0.5*DL, 0.5*DW, 0.5*DH);
		// Remove parts at the right part of the cliff
		Vecd halfsize_linner(0.5*VWx, 0.5*DW, 0.5*DH);
		Vecd translation_linner(0.5*VWx, 0.5*DW, 0.5*DH);
		// Remove parts above cliff
		Vecd halfsize_rinner(0.5*(DL - VWx), 0.5*DW, 0.5*(DH - VWH));
		Vecd translation_rinner(VWx + 0.5*(DL - VWx), 0.5*DW, VWH + 0.5*(DH - VWH));

		body_shape_.add<TriangleMeshShapeBrick>(halfsize_outer, resolution, translation_outer);
		body_shape_.substract<TriangleMeshShapeBrick>(halfsize_linner, resolution, translation_linner);
		body_shape_.substract<TriangleMeshShapeBrick>(translation_rinner, resolution, translation_rinner);

		// remove inside of cliff
		Vecd halfsize_wall_inner(0.5*(DL - VWx - BW), 0.5*DW + BW, 0.5*VWH);
		Vecd translation_wall_inner(0.5*(DL + VWx + BW), 0.5*DW, 0.5*VWH - BW);
		body_shape_.substract<TriangleMeshShapeBrick>(halfsize_wall_inner, resolution, translation_wall_inner);
	}
};

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

class Boulder : public SolidBody
{
public:
	Boulder(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		Vecd halfsize(0.5 * BDL, 0.5 * BDW, 0.5 * BDH);
		Vecd translation_wall(VWx - BDx - 0.5*BDL, BDy, BDz + 0.5*BDH);

		body_shape_.add<TriangleMeshShapeBrick>(halfsize, 0, translation_wall);
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