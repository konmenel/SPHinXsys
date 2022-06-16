#pragma once
#include "sphinxsys.h"

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemNSC.h"

namespace SPH 
{
	using chrono::ChSystem;
	using chrono::ChSystemNSC;
	using chrono::ChVector;
	using chrono::Quaternion;
	using chrono::ChFrameMoving;
	using chrono::ChBody;
	using chrono::ChBodyEasyBox;
	using chrono::ChMaterialSurfaceNSC;
	using chrono::ChForce;
	using chrono::ChTimestepper;
	using chrono::ChSolver;

	// Converts Vec3d to chrono::ChVector<>
	inline ChVector<> vecToCh(const Vec3d &vec) { return ChVector<>(vec[0], vec[1], vec[2]); }

	// Converts chrono::ChVector<> to Vec3d
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
				initial_frame_loc_ = vecToSim(ch_body->GetPos());
			}
			
		virtual ~ConstrainSolidBodyPartByChrono(){};

	protected:
		std::shared_ptr<ChBody> ch_body_;
		Vec3d initial_frame_loc_;

		virtual void setupDynamics(Real dt = 0.0) override 
		{
			body_->setNewlyUpdated();
		}

		virtual void Update(size_t index_i, Real dt = 0.0) override 
		{
			ChVector<> rpos, pos, vel, acc;
			Quaternion rot;
			
			// Can calculate this once for all particle in constuctor and
			// store as attribute. Increase speed in exchange for memory.
			rpos = vecToCh(pos_0_[index_i] - initial_frame_loc_);	// Position with respect to local frame

			pos = ch_body_->TransformDirectionLocalToParent(rpos); 	// Translation of local frame
			vel = ch_body_->PointSpeedLocalToParent(rpos); 			// Velocity of local frame
			acc = ch_body_->PointAccelerationLocalToParent(rpos); 	// Velocity of local frame
			rot = ch_body_->GetRot();

			pos_n_[index_i] = vecToSim(pos);
			vel_n_[index_i] = vecToSim(vel);
			dvel_dt_[index_i] = vecToSim(acc);
			n_[index_i] = vecToSim(
				rot.Rotate(
					vecToCh(n_0_[index_i])
				)
			);
		}
	};


	/**
	 * @class TotalForceOnSolidBodyPartForChrono
	 * @brief Compute the force acting on the solid body part
	 * for applying to chrono forces latter
	 */
	class TotalForceOnSolidBodyPartForChrono
		: public PartDynamicsByParticleReduce<SimTK::SpatialVec, ReduceSum<SimTK::SpatialVec>>,
		public solid_dynamics::SolidDataSimple
	{
	public:
		TotalForceOnSolidBodyPartForChrono(SolidBody &solid_body,
										BodyRegionByParticle &body_part,
										std::shared_ptr<ChBody> ch_body,
										ChSystem &ch_system)
		: PartDynamicsByParticleReduce<SimTK::SpatialVec, ReduceSum<SimTK::SpatialVec>>(solid_body, body_part),
		solid_dynamics::SolidDataSimple(solid_body),
		force_from_fluid_(particles_->force_from_fluid_), contact_force_(particles_->contact_force_),
		pos_n_(particles_->pos_n_), ch_system_(ch_system), ch_body_(ch_body) 
		{
			initial_reference_ = SimTK::SpatialVec(Vec3d(0), Vec3d(0));
		}
		
		virtual ~TotalForceOnSolidBodyPartForChrono(){};

	protected:
		StdLargeVec<Vecd> &force_from_fluid_, &contact_force_, &pos_n_;
		ChSystem &ch_system_;
		std::shared_ptr<ChBody> ch_body_;
		ChVector<> current_mobod_origin_location_;

		virtual void SetupReduce() override
		{
			current_mobod_origin_location_ = ch_body_->GetPos();
		}

		virtual SimTK::SpatialVec ReduceFunction(size_t index_i, Real dt = 0.0) override
		{
			// TODO: Maybe try accumulate total force each timestep (remember to clear though)
			// Calculate force
			Vec3d force_from_particle = force_from_fluid_[index_i] + contact_force_[index_i];
			
			// Calculate torque
			Vec3d displacement = pos_n_[index_i] - vecToSim(current_mobod_origin_location_);
			Vec3d torque_from_particle = cross(displacement, force_from_particle);

			return SimTK::SpatialVec(torque_from_particle, force_from_particle);
		}
	};
} // namespace SPH
