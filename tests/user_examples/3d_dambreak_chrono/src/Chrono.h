#include "sphinxsys.h"

#define SMOOTH_CONTACT 0

#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChBodyEasy.h"

#if SMOOTH_CONTACT
#include "chrono/physics/ChSystemSMC.h"
using SystemCh = chrono::ChSystemSMC;
using SurfMaterialCh = chrono::ChMaterialSurfaceSMC;
#else
#include "chrono/physics/ChSystemNSC.h"
using SystemCh = chrono::ChSystemNSC;
using SurfMaterialCh = chrono::ChMaterialSurfaceNSC;
#endif // end if SMOOTH_CONTACT

using namespace SPH;

using chrono::ChSystem;
using chrono::ChVector;
using chrono::Quaternion;
using chrono::ChBody;
using chrono::ChBodyEasyBox;
using chrono::ChTimestepper;
using chrono::ChSolver;

const Real resolution_ref = 0.01;	  	//particle spacing
const Real DL = 1.0;
const Real DW = 1.0;
const Real DH = 1.0;
const Real BDL = 0.1;
const Real BDW = 0.1;
const Real BDH = 0.1;
const Real Bx = DL * 0.5;
const Real By = DW * 0.5;
const Real Bz = DH * 0.5;
const Real BW = 4 * resolution_ref;

const BoundingBox system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DL + BW, DW + BW, DH + BW));


const Real rho0_s = 2.8e3;								/**< Boulder Density [kg/m^3] from paper. */
const Real boulder_vol = BDL * BDH * BDW;				/**< Boulder Volume [m^3]. (1.5 x 2.0 x 3.0 cm) */
const Real boulder_mass = rho0_s * boulder_vol;			/**< Boulder Mass [kg]. */
const Real poisson = 0.3;								/**< Poisson's ratio. */
const Real Youngs_modulus = 73e9;						/**< Young's modulus [Pa]. */
const Real surface_thickness = 1.0;

const int resolution = 20;

auto collition_model = chrono_types::make_shared<chrono::collision::ChCollisionModelBullet>();

inline ChVector<> vecToCh(const Vec3d &vec) { return ChVector<>(vec[0], vec[1], vec[2]); }
inline Vec3d vecToSim(const ChVector<> &vec) { return Vec3d(vec[0], vec[1], vec[2]); }

std::shared_ptr<ChBody> addBoxCh(ChSystem &ch_system)
{
	auto box_mat = chrono_types::make_shared<SurfMaterialCh>();
	box_mat->SetFriction(0.4f);
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
	box_ch->SetPos(ChVector<>(Bx, By, Bz));

	ch_system.AddBody(box_ch);
	return box_ch;
}

void addWallsCh(ChSystem &ch_system)
{	
	auto box_mat = chrono_types::make_shared<SurfMaterialCh>();
	box_mat->SetFriction(0.2f);
#if SMOOTH_CONTACT
	box_mat->SetYoungModulus(Youngs_modulus);
	box_mat->SetPoissonRatio(poisson);
#endif // end if SMOOTH_CONTACT

	auto floor_ch = chrono_types::make_shared<ChBodyEasyBox>(	DL,
																DW,
																BW,
																rho0_s,
																false,
																true,
																box_mat);
	floor_ch->SetName("Floor");
	floor_ch->SetPos(ChVector<>(DL * 0.5, DW * 0.5, -BW * 0.5));
	floor_ch->SetBodyFixed(true);

	ch_system.AddBody(floor_ch);
	return;
}

class WallBoundary : public SolidBody
{
public:
	WallBoundary(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		Vecd halfsize_outer(0.5*DL + BW, 0.5*DW + BW, 0.5*DH + BW);
		Vecd translation_outer(0.5*DL, 0.5*DW, 0.5*DH);
		Vecd halfsize_inner(0.5*DL + BW, 0.5*DW, 0.5*DH);
		Vecd translation_inner(0.5*DL, 0.5*DW, 0.5*DH);

		body_shape_.add<TriangleMeshShapeBrick>(halfsize_outer, resolution, translation_outer);
		body_shape_.substract<TriangleMeshShapeBrick>(halfsize_inner, resolution, translation_inner);
	}
};

class Box : public SolidBody
{
public:
	Box(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		Vecd halfsize(0.5*BDL, 0.5*BDW, 0.5*BDH);
		Vecd translation(Bx, By, Bz);

		body_shape_.add<TriangleMeshShapeBrick>(halfsize, resolution, translation);
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

// TODO: Implement this class
class ConstrainSolidBodyPartByChrono : public solid_dynamics::ConstrainSolidBodyRegion
{
public:
	ConstrainSolidBodyPartByChrono(	SolidBody &solid_body,
									BodyRegionByParticle &body_part,
									std::shared_ptr<ChBody> ch_body)
		: ConstrainSolidBodyRegion(solid_body, body_part), ch_body_(ch_body)
		{
			initial_body_origin_location_ = vecToSim(ch_body->GetPos());
		}
		
	virtual ~ConstrainSolidBodyPartByChrono(){};

protected:
	std::shared_ptr<ChBody> ch_body_;
	Vec3d initial_body_origin_location_;

	virtual void setupDynamics(Real dt = 0.0) override 
	{
		body_->setNewlyUpdated();
	}

	void virtual Update(size_t index_i, Real dt = 0.0) override 
	{
		Vec3d rr;
		ChVector<> pos, vel, acc;
		Quaternion rot;
		rr = pos_0_[index_i] - initial_body_origin_location_;

		ChVector<> rr_ch = vecToCh(rr);
		
		pos = ch_body_->GetPos(); 							// Translation of origin
		vel = ch_body_->PointSpeedLocalToParent(rr_ch); 		// Velocity of origin
		acc = ch_body_->PointAccelerationLocalToParent(rr_ch); 	// Velocity of point
		rot = ch_body_->GetRot();								// Rotation of local fram

		pos_n_[index_i] = vecToSim(pos) + rr;
		vel_n_[index_i] = vecToSim(vel);
		dvel_dt_[index_i] = vecToSim(acc);
		n_[index_i] = vecToSim(rot.Rotate(vecToCh(n_0_[index_i])));
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