#pragma once

// Removing MSVC compiler min/max macros because of clash
#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

#include "solid_dynamics.h"

#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemNSC.h"

#include "H5File.hpp"

using chrono::ChSystemNSC;
using chrono::ChVector;
using chrono::Quaternion;
using chrono::ChBody;
using chrono::ChBodyEasyBox;
using chrono::ChMaterialSurfaceNSC;

namespace SPH 
{
	// Converts Vec3d to chrono::ChVector<>
	inline ChVector<> vecToCh(const Vec3d &vec) { return ChVector<>(vec[0], vec[1], vec[2]); }

	// Converts chrono::ChVector<> to Vec3d
	inline Vec3d vecToSim(const ChVector<> &vec) { return Vec3d(vec[0], vec[1], vec[2]); }


	class SolidBodyWithChrono : public SolidBody
	{
	public:
		SolidBodyWithChrono(SPHSystem &system, const std::string &body_name,
							ChSystemNSC &ch_system)
			: SolidBody(system, body_name), ch_system_(ch_system) {}

		virtual ~SolidBodyWithChrono(){}

		virtual void setChronoBody(std::shared_ptr<ChBody> ch_body) 
		{
			ch_system_.AddBody(ch_body);
			ch_body_ = ch_body;
		}

		virtual ChSystemNSC& getChSystem() { return ch_system_; }

		virtual std::shared_ptr<ChBody> getChBody() { return ch_body_; }
		
		virtual void writeChronoDataToH5(HF::Group &parent) 
		{
			const ChVector<> &pos = ch_body_->GetPos();
			const ChVector<> &vel = ch_body_->GetPos_dt();
			const ChVector<> &acc = ch_body_->GetPos_dtdt();

			const ChVector<> angles = ch_body_->GetRot().Q_to_Euler123();	// ChVector(roll, pitch, yaw)
			const ChVector<> angvel = ch_body_->GetWvel_par();
			const ChVector<> angacc = ch_body_->GetWacc_par();
			
			std::vector<size_t> dims = { 3 };
			HF::DataSpace dims_dataspace(dims);

			HF::DataSet pos_dataset = parent.createDataSet<Real>("Pos", dims_dataspace);
			pos_dataset.write(&pos[0]);

			HF::DataSet vel_dataset = parent.createDataSet<Real>("Vel", dims_dataspace);
			vel_dataset.write(&vel[0]);
			
			HF::DataSet acc_dataset = parent.createDataSet<Real>("Acc", dims_dataspace);
			acc_dataset.write(&acc[0]);
			
			HF::DataSet ang_dataset = parent.createDataSet<Real>("Angles", dims_dataspace);
			ang_dataset.write(&angles[0]);
			
			HF::DataSet angvel_dataset = parent.createDataSet<Real>("AngVel", dims_dataspace);
			angvel_dataset.write(&angvel[0]);
			
			HF::DataSet angacc_dataset = parent.createDataSet<Real>("AngAcc", dims_dataspace);
			angacc_dataset.write(&angacc[0]);
		}
	protected:
		ChSystemNSC &ch_system_;
		std::shared_ptr<ChBody> ch_body_;
	};


	namespace solid_dynamics
	{
		typedef DataDelegateSimple<SolidBodyWithChrono, SolidParticles, Solid> SolidChronoDataSimple;

		class ConstrainSolidChronoBodyRegion : public PartSimpleDynamicsByParticle, public SolidChronoDataSimple
		{
		public:
			ConstrainSolidChronoBodyRegion(SPHBody &sph_body, BodyPartByParticle &body_part)
			: PartSimpleDynamicsByParticle(sph_body, body_part), SolidChronoDataSimple(sph_body),
			  pos_n_(particles_->pos_n_), pos_0_(particles_->pos_0_),
			  n_(particles_->n_), n_0_(particles_->n_0_),
			  vel_n_(particles_->vel_n_), dvel_dt_(particles_->dvel_dt_),
			  vel_ave_(particles_->vel_ave_), dvel_dt_ave_(particles_->dvel_dt_ave_) {}
			virtual ~ConstrainSolidChronoBodyRegion(){}

		protected:
			StdLargeVec<Vecd> &pos_n_, &pos_0_;
			StdLargeVec<Vecd> &n_, &n_0_;
			StdLargeVec<Vecd> &vel_n_, &dvel_dt_, &vel_ave_, &dvel_dt_ave_;
			virtual Vecd getDisplacement(Vecd &pos_0, Vecd &pos_n) { return pos_n; };
			virtual Vecd getVelocity(Vecd &pos_0, Vecd &pos_n, Vecd &vel_n) { return Vecd(0); };
			virtual Vecd getAcceleration(Vecd &pos_0, Vecd &pos_n, Vecd &dvel_dt) { return Vecd(0); };
			virtual SimTK::Rotation getBodyRotation(Vecd &pos_0, Vecd &pos_n, Vecd &dvel_dt) { return SimTK::Rotation(); }
			virtual void Update(size_t index_i, Real dt = 0.0) override
			{
				Vecd pos_0 = pos_0_[index_i];
				Vecd pos_n = pos_n_[index_i];
				Vecd vel_n = vel_n_[index_i];
				Vecd dvel_dt = dvel_dt_[index_i];

				pos_n_[index_i] = getDisplacement(pos_0, pos_n);
				vel_n_[index_i] = getVelocity(pos_0, pos_n, vel_n);
				dvel_dt_[index_i] = getAcceleration(pos_0, pos_n, dvel_dt);
				/** the average values are prescirbed also. */
				vel_ave_[index_i] = vel_n_[index_i];
				dvel_dt_ave_[index_i] = dvel_dt_[index_i];
			}
		};
	}


	/**
	 * @class ConstrainSolidBodyPartByChrono
	 * @brief Constrain a solid body part from the motion
	 * computed from Chrono.
	 */
	class ConstrainSolidBodyPartByChrono : public solid_dynamics::ConstrainSolidChronoBodyRegion
	{
	public:
		ConstrainSolidBodyPartByChrono(	SolidBodyWithChrono &solid_body,
										BodyRegionByParticle &body_part)
			: ConstrainSolidChronoBodyRegion(solid_body, body_part), ch_body_(body_->getChBody())
			{
				initial_frame_loc_ = vecToSim(ch_body_->GetPos());
			}
			
		virtual ~ConstrainSolidBodyPartByChrono(){};

	protected:
		Vec3d initial_frame_loc_;
		std::shared_ptr<ChBody> ch_body_;

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

			pos = ch_body_->GetPos();
			vel = ch_body_->PointSpeedLocalToParent(rpos);
			acc = ch_body_->PointAccelerationLocalToParent(rpos);
			rot = ch_body_->GetRot();

			pos_n_[index_i] = vecToSim(pos + rot.Rotate(rpos));
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
	 * for applying to chrono forces latter.
	 */
	class TotalForceOnSolidBodyPartForChrono
		: public PartDynamicsByParticleReduce<SimTK::SpatialVec, ReduceSum<SimTK::SpatialVec>>,
		public solid_dynamics::SolidChronoDataSimple
	{
	public:
		TotalForceOnSolidBodyPartForChrono(SolidBodyWithChrono &solid_body,
										BodyRegionByParticle &body_part)
		: PartDynamicsByParticleReduce<SimTK::SpatialVec, ReduceSum<SimTK::SpatialVec>>(solid_body, body_part),
		solid_dynamics::SolidChronoDataSimple(solid_body),
		force_from_fluid_(particles_->force_from_fluid_), contact_force_(particles_->contact_force_),
		pos_n_(particles_->pos_n_), ch_system_(body_->getChSystem()), ch_body_(body_->getChBody()) 
		{
			initial_reference_ = SimTK::SpatialVec(Vec3d(0), Vec3d(0));
		}
		
		virtual ~TotalForceOnSolidBodyPartForChrono(){};

	protected:
		StdLargeVec<Vec3d> &force_from_fluid_, &contact_force_, &pos_n_;
		ChSystemNSC &ch_system_;
		std::shared_ptr<ChBody> ch_body_;
		Vec3d current_mobod_origin_location_;

		virtual void SetupReduce() override
		{
			current_mobod_origin_location_ = vecToSim(ch_body_->GetPos());
		}

		virtual SimTK::SpatialVec ReduceFunction(size_t index_i, Real dt = 0.0) override
		{
			// TODO: Maybe try accumulate total force each timestep (remember to clear though)
			// Calculate force
			Vec3d force_from_particle = force_from_fluid_[index_i] + contact_force_[index_i];
			
			// Calculate torque
			Vec3d displacement = pos_n_[index_i] - current_mobod_origin_location_;
			Vec3d torque_from_particle = cross(displacement, force_from_particle);

			return SimTK::SpatialVec(torque_from_particle, force_from_particle);
		}
	};

	/**
	 * @class WriteChronoBodyData
	 * @brief Writes the Chrono data to a csv file.
	 */
	class WriteChronoBodyData
	{
	public:
		WriteChronoBodyData(In_Output &in_output,
							ChSystemNSC &ch_system,
							std::shared_ptr<ChBody> ch_body)
		: in_output_(in_output), ch_system_(ch_system), ch_body_(ch_body),
		  filefullpath_(in_output_.output_folder_ + "/chrono_" + ch_body_->GetName() + "_data.csv")
		{
			std::ofstream out_file(filefullpath_.c_str(), std::ios::app);
			// header
			out_file << "#THE UNITS PRESENTED BELOW ARE JUST FOR REFERENCE. THE TRUE UNITS ARE " 
					<< "DEPENDENT ON THE UNITS USED TO SET UP THE SIMULATION\n";
			
			out_file << "time [s],"
					<< "COM.x [m],"
					<< "COM.y [m],"
					<< "COM.z [m],"
					<< "vel.x [m/s],"
					<< "vel.y [m/s],"
					<< "vel.z [m/s],"
					<< "acc.x [m/s^2],"
					<< "acc.y [m/s^2],"
					<< "acc.z [m/s^2],"
					<< "rot.roll [rad],"
					<< "rot.pitch [rad],"
					<< "rot.yaw [rad],"
					<< "angvel.x [rad/s],"
					<< "angvel.y [rad/s],"
					<< "angvel.z [rad/s],"
					<< "angacc.x [rad/s^2],"
					<< "angacc.y [rad/s^2],"
					<< "angacc.z [rad/s^2]"
					<< "\n";

			out_file.close();
		}

		virtual ~WriteChronoBodyData(){}

		virtual void writeToFile(const size_t iteration_step = 0)
		{
			std::ofstream out_file(filefullpath_.c_str(), std::ios::app);
			
			const ChVector<> &pos = ch_body_->GetPos();
			const ChVector<> &vel = ch_body_->GetPos_dt();
			const ChVector<> &acc = ch_body_->GetPos_dtdt();

			const ChVector<> angles = ch_body_->GetRot().Q_to_Euler123();	// ChVector(roll, pitch, yaw)
			const ChVector<> angvel = ch_body_->GetWvel_par();
			const ChVector<> angacc = ch_body_->GetWacc_par();

			out_file << GlobalStaticVariables::physical_time_ << ","

					<< pos[0] << ","
					<< pos[1] << ","
					<< pos[2] << ","

					<< vel[0] << ","
					<< vel[1] << ","
					<< vel[2] << ","

					<< acc[0] << ","
					<< acc[1] << ","
					<< acc[2] << ","

					<< angles[0] << ","
					<< angles[1] << ","
					<< angles[2] << ","

					<< angvel[0] << ","
					<< angvel[1] << ","
					<< angvel[2] << ","

					<< angacc[0] << ","
					<< angacc[1] << ","
					<< angacc[2] << "\n";

			out_file.close();
		}
	
	protected:
		std::string filefullpath_;
		In_Output &in_output_;
		ChSystemNSC &ch_system_;
		std::shared_ptr<ChBody> ch_body_;
	};
} // namespace SPH
