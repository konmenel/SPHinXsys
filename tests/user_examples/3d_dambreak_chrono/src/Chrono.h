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
using chrono::ChForce;
using chrono::ChTimestepper;
using chrono::ChSolver;

// TODO: Add the rest of the walls to chronos
// TODO: True to create each wall of the tank seperately
// TODO: Change tank materials
// TODO: Tune contact coefficients


inline ChVector<> vecToCh(const Vec3d &vec) { return ChVector<>(vec[0], vec[1], vec[2]); }
inline Vec3d vecToSim(const ChVector<> &vec) { return Vec3d(vec[0], vec[1], vec[2]); }

/**
 * @class ConstrainSolidBodyPartByChrono
 * @brief Constrain a solid body part from the motion
 * computed from Chrono.
 */
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
 * @class TotalForceOnSolidBodyPartForChrono
 * @brief Compute the force acting on the solid body part
 * for applying to chrono forces latter
 */
class TotalForceOnSolidBodyPartForChrono
	: public PartSimpleDynamicsByParticle, public solid_dynamics::SolidDataSimple
{
public:
	TotalForceOnSolidBodyPartForChrono(SolidBody &solid_body,
										BodyRegionByParticle &body_part,
										ChSystem &ch_system,
										std::shared_ptr<ChBody> ch_body)
	: PartSimpleDynamicsByParticle(solid_body, body_part),
	  solid_dynamics::SolidDataSimple(solid_body),
	  force_from_fluid_(particles_->force_from_fluid_), contact_force_(particles_->contact_force_),
	  pos_n_(particles_->pos_n_), ch_system_(ch_system), ch_body_(ch_body) {}
	
	virtual ~TotalForceOnSolidBodyPartForChrono(){};

protected:
	StdLargeVec<Vecd> &force_from_fluid_, &contact_force_, &pos_n_;
	ChSystem &ch_system_;
	std::shared_ptr<ChBody> ch_body_;
	ChVector<> current_mobod_origin_location_;

	virtual void setupDynamics(Real dt = 0.0) override
	{
		current_mobod_origin_location_ = ch_body_->GetPos();
		ch_body_->RemoveAllForces();
	}

	virtual void Update(size_t index_i, Real dt = 0.0) override
	{
		// Calculate force
		Vecd force_from_particle = force_from_fluid_[index_i] + contact_force_[index_i];
		
		// Calculate torque
		Vecd displacement = pos_n_[index_i] - vecToSim(current_mobod_origin_location_);
		Vecd torque_from_particle = cross(displacement, force_from_particle);

		ch_body_->Accumulate_force(vecToCh(force_from_particle),
								   ChVector<>(0, 0, 0),
								   false);
		ch_body_->Accumulate_torque(vecToCh(torque_from_particle),
									false);
	}
};
