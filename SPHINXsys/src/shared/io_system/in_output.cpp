/**
 * @file 	in_output.cpp
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#include "in_output.h"
#include "all_bodies.h"
#include "level_set.h"
#include "sph_system.h"

namespace SPH
{
	//=============================================================================================//
	In_Output::In_Output(SPHSystem &sph_system)
		: sph_system_(sph_system),
		  input_folder_("./input"), output_folder_("./output"),
		  restart_folder_("./restart"), reload_folder_("./reload")
	{
		if (!fs::exists(input_folder_))
		{
			fs::create_directory(input_folder_);
		}

		if (!fs::exists(output_folder_))
		{
			fs::create_directory(output_folder_);
		}

		if (!fs::exists(restart_folder_))
		{
			fs::create_directory(restart_folder_);
		}

		if (sph_system.restart_step_ == 0)
		{
			fs::remove_all(restart_folder_);
			fs::create_directory(restart_folder_);

			fs::remove_all(output_folder_);
			fs::create_directory(output_folder_);
		}

		restart_step_ = std::to_string(sph_system.restart_step_);

		sph_system.in_output_ = this;
	}
	//=============================================================================================//
	void PltEngine::
		writeAQuantityHeader(std::ofstream &out_file, const Real &quantity, const std::string &quantity_name)
	{
		out_file << "\"" << quantity_name << "\""
				 << "   ";
	}
	//=============================================================================================//
	void PltEngine::
		writeAQuantityHeader(std::ofstream &out_file, const Vecd &quantity, const std::string &quantity_name)
	{
		for (int i = 0; i != Dimensions; ++i)
			out_file << "\"" << quantity_name << "[" << i << "]\""
					 << "   ";
	}
	//=============================================================================================//
	void PltEngine::writeAQuantity(std::ofstream &out_file, const Real &quantity)
	{
		out_file << std::fixed << std::setprecision(9) << quantity << "   ";
	}
	//=============================================================================================//
	void PltEngine::writeAQuantity(std::ofstream &out_file, const Vecd &quantity)
	{
		for (int i = 0; i < Dimensions; ++i)
			out_file << std::fixed << std::setprecision(9) << quantity[i] << "   ";
	}
	//=============================================================================================//
	std::string BodyStatesIO::convertPhysicalTimeToString(Real convertPhysicalTimeToStream)
	{
		int i_time = int(GlobalStaticVariables::physical_time_ * 1.0e6);
		std::stringstream s_time;
		s_time << std::setw(10) << std::setfill('0') << i_time;
		return s_time.str();
	}
	//=============================================================================================//
	void BodyStatesRecordingToVtp::writeWithFileName(const std::string &sequence)
	{
		for (SPHBody *body : bodies_)
		{
			if (body->checkNewlyUpdated())
			{
				//TODO: we can short the file name by without using SPHBody
				std::string filefullpath = in_output_.output_folder_ + "/SPHBody_" + body->getBodyName() + "_" + sequence + ".vtp";
				if (fs::exists(filefullpath))
				{
					fs::remove(filefullpath);
				}
				std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);
				//begin of the XML file
				out_file << "<?xml version=\"1.0\"?>\n";
				out_file << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
				out_file << " <PolyData>\n";

				BaseParticles *base_particles = body->base_particles_;
				size_t total_real_particles = base_particles->total_real_particles_;
				out_file << "  <Piece Name =\"" << body->getBodyName() << "\" NumberOfPoints=\"" << total_real_particles << "\" NumberOfVerts=\"" << total_real_particles << "\">\n";

				body->writeParticlesToVtpFile(out_file);

				out_file << "   </PointData>\n";

				//write empty cells
				out_file << "   <Verts>\n";
				out_file << "    <DataArray type=\"Int32\"  Name=\"connectivity\"  Format=\"ascii\">\n";
				out_file << "    ";
				for (size_t i = 0; i != total_real_particles; ++i)
				{
					out_file << i << " ";
				}
				out_file << std::endl;
				out_file << "    </DataArray>\n";
				out_file << "    <DataArray type=\"Int32\"  Name=\"offsets\"  Format=\"ascii\">\n";
				out_file << "    ";
				for (size_t i = 0; i != total_real_particles; ++i)
				{
					out_file << i + 1 << " ";
				}
				out_file << std::endl;
				out_file << "    </DataArray>\n";
				out_file << "   </Verts>\n";

				out_file << "  </Piece>\n";

				out_file << " </PolyData>\n";
				out_file << "</VTKFile>\n";

				out_file.close();
			}
			body->setNotNewlyUpdated();
		}
	}
	//=============================================================================================//
	void BodyStatesRecordingToLegacyVtk::writeWithFileName(const std::string &sequence)
	{
		for (SPHBody *body : bodies_)
		{
			if (body->checkNewlyUpdated())
			{
				//TODO: we can short the file name by without using SPHBody
				std::string filefullpath = in_output_.output_folder_ + "/SPHBody_" + body->getBodyName() + "_" + sequence + ".vtk";
				if (fs::exists(filefullpath))
				{
					fs::remove(filefullpath);
				}
				std::ofstream out_file(filefullpath.c_str(), std::ios::trunc | std::ios::binary);
				//begin of the Legacy VTK file
				out_file << "# vtk DataFile Version 3.0\n";
				out_file << body->getBodyName() << " output\n";
				out_file << "BINARY\n";
				out_file << "DATASET POLYDATA\n";

				BaseParticles *base_particles = body->base_particles_;
				size_t total_real_particles = base_particles->total_real_particles_;
				
				//write the partcle data
				body->writeParticlesToBinLecVtkFile(out_file);
				
				out_file.close();
			}
			body->setNotNewlyUpdated();
		}
	}
	//=============================================================================================//
	void BodyStatesRecordingToVtuString::writeWithFileName(const std::string& sequence)
	{
		for (SPHBody* body : bodies_)
		{
			if (body->checkNewlyUpdated())
			{
				const auto& vtuName = body->getBodyName() + "_" + sequence + ".vtu";
				std::stringstream sstream;
				//begin of the XML file
				writeVtu(sstream, body);
				_vtuData[vtuName] = sstream.str();
			}
			body->setNotNewlyUpdated();
		}
	}
	//=============================================================================================//
	void BodyStatesRecordingToVtuString::writeVtu(std::ostream& stream, SPHBody* body) const
	{
		stream << "<?xml version=\"1.0\"?>\n";
		stream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
		stream << " <UnstructuredGrid>\n";

		BaseParticles* base_particles = body->base_particles_;
		size_t total_real_particles = base_particles->total_real_particles_;
		stream << "  <Piece Name =\"" << body->getBodyName() << "\" NumberOfPoints=\"" << total_real_particles << "\" NumberOfCells=\"0\">\n";

		body->writeParticlesToVtuFile(stream);

		stream << "   </PointData>\n";

		//write empty cells
		stream << "   <Cells>\n";
		stream << "    <DataArray type=\"Int32\"  Name=\"connectivity\"  Format=\"ascii\">\n";
		stream << "    </DataArray>\n";
		stream << "    <DataArray type=\"Int32\"  Name=\"offsets\"  Format=\"ascii\">\n";
		stream << "    </DataArray>\n";
		stream << "    <DataArray type=\"types\"  Name=\"offsets\"  Format=\"ascii\">\n";
		stream << "    </DataArray>\n";
		stream << "   </Cells>\n";

		stream << "  </Piece>\n";

		stream << " </UnstructuredGrid>\n";
		stream << "</VTKFile>\n";
	}
	//=============================================================================================//
	const BodyStatesRecordingToVtuString::VtuStringData& BodyStatesRecordingToVtuString::GetVtuData() const
	{
		return _vtuData;
	}
	//=============================================================================================//
	SurfaceOnlyBodyStatesRecordingToVtu::SurfaceOnlyBodyStatesRecordingToVtu(In_Output& in_output, SPHBodyVector bodies)
			: BodyStatesRecording(in_output, bodies),
			surface_body_layer_vector_({})
	{
		for (SPHBody* body : bodies_) surface_body_layer_vector_.push_back(BodySurface(*body));
	}
	//=============================================================================================//
	void SurfaceOnlyBodyStatesRecordingToVtu::writeWithFileName(const std::string& sequence)
	{
		for (size_t i = 0; i < bodies_.size(); i++)
		//for (size_t i = 0; surface_body_layer_vector_.size(); i++)
		{
			SPHBody* body = bodies_[i];
			if (body->checkNewlyUpdated())
			{
				std::string filefullpath = in_output_.output_folder_ + "/SPHBody_" + body->getBodyName() + "_" + sequence + ".vtu";
				if (fs::exists(filefullpath))
				{
					fs::remove(filefullpath);
				}
				std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);
				//begin of the XML file
				out_file << "<?xml version=\"1.0\"?>\n";
				out_file << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
				out_file << " <PolyData>\n";

				BaseParticles *base_particles = body->base_particles_;
				size_t total_real_particles = base_particles->total_real_particles_;
				out_file << "  <Piece Name =\"" << body->getBodyName() << "\" NumberOfPoints=\"" << total_real_particles << "\" NumberOfVerts=\"" << total_real_particles << "\">\n";

				body->writeParticlesToVtpFile(out_file);

				out_file << "   </PointData>\n";

				//write empty cells
				out_file << "   <Verts>\n";
				out_file << "    <DataArray type=\"Int32\"  Name=\"connectivity\"  Format=\"ascii\">\n";
				out_file << "    ";
				for (size_t i = 0; i != total_real_particles; ++i)
				{
					out_file << i << " ";
				}
				out_file << std::endl;
				out_file << "    </DataArray>\n";
				out_file << "    <DataArray type=\"Int32\"  Name=\"offsets\"  Format=\"ascii\">\n";
				out_file << "    ";
				for (size_t i = 0; i != total_real_particles; ++i)
				{
					out_file << i + 1 << " ";
				}
				out_file << std::endl;
				out_file << "    </DataArray>\n";
				out_file << "   </Verts>\n";

				out_file << "  </Piece>\n";

				out_file << " </PolyData>\n";
				out_file << "</VTKFile>\n";

				out_file.close();
			}
			body->setNotNewlyUpdated();
		}
	}
	//=============================================================================================//
	void BodyStatesRecordingToPlt::writeWithFileName(const std::string &sequence)
	{
		for (SPHBody *body : bodies_)
		{
			if (body->checkNewlyUpdated())
			{
				std::string filefullpath = in_output_.output_folder_ + "/SPHBody_" + body->getBodyName() + "_" + sequence + ".plt";
				if (fs::exists(filefullpath))
				{
					fs::remove(filefullpath);
				}
				std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);

				//begin of the plt file writing

				body->writeParticlesToPltFile(out_file);

				out_file.close();
			}
			body->setNotNewlyUpdated();
		}
	}
	//=============================================================================================//
	WriteToVtpIfVelocityOutOfBound::
		WriteToVtpIfVelocityOutOfBound(In_Output &in_output, SPHBodyVector bodies, Real velocity_bound)
		: BodyStatesRecordingToVtp(in_output, bodies), out_of_bound_(false)
	{
		std::transform(bodies.begin(), bodies.end(), std::back_inserter(check_bodies_),
					   [&](SPHBody *body) -> VelocityBoundCheck
					   { return VelocityBoundCheck(*body, velocity_bound); });
	}
	//=============================================================================================//
	void WriteToVtpIfVelocityOutOfBound::writeWithFileName(const std::string &sequence)
	{
		for (auto check_body : check_bodies_)
		{
			out_of_bound_ = out_of_bound_ || check_body.parallel_exec();
		}

		if (out_of_bound_)
		{
			BodyStatesRecordingToVtp::writeWithFileName(sequence);
			std::cout << "\n Velocity is out of bound at iteration step " << sequence
					  << "\n The body states have been outputted and the simulation terminates here. \n";
		}
	}
	//=============================================================================================//
	MeshRecordingToPlt ::MeshRecordingToPlt(In_Output &in_output, SPHBody &body, BaseMeshField *mesh_field)
		: BodyStatesRecording(in_output, body), mesh_field_(mesh_field)
	{
		filefullpath_ = in_output_.output_folder_ + "/" + body.getBodyName() + "_" + mesh_field_->Name() + ".dat";
	}
	//=============================================================================================//
	void MeshRecordingToPlt::writeWithFileName(const std::string &sequence)
	{
		std::ofstream out_file(filefullpath_.c_str(), std::ios::app);
		mesh_field_->writeMeshFieldToPlt(out_file);
		out_file.close();
	}
	//=============================================================================================//
	ReloadParticleIO::ReloadParticleIO(In_Output &in_output, SPHBodyVector bodies) : BodyStatesIO(in_output, bodies)
	{
		if (!fs::exists(in_output.reload_folder_))
		{
			fs::create_directory(in_output.reload_folder_);
		}

		std::transform(bodies.begin(), bodies.end(), std::back_inserter(file_paths_),
					   [&](SPHBody *body) -> std::string
					   { return in_output.reload_folder_ + "/SPHBody_" + body->getBodyName() + "_rld.xml"; });
	};
	//=============================================================================================//
	ReloadParticleIO::ReloadParticleIO(In_Output &in_output, SPHBodyVector bodies,
									   StdVec<std::string> given_body_names) : ReloadParticleIO(in_output, bodies)
	{
		std::transform(given_body_names.begin(), given_body_names.end(), file_paths_.begin(),
					   [&](const std::string &body_name) -> std::string
					   { return in_output.reload_folder_ + "/SPHBody_" + body_name + "_rld.xml"; });
	}
	//=============================================================================================//
	void ReloadParticleIO::writeToFile(size_t iteration_step)
	{
		for (size_t i = 0; i < bodies_.size(); ++i)
		{
			std::string filefullpath = file_paths_[i];

			if (fs::exists(filefullpath))
			{
				fs::remove(filefullpath);
			}
			bodies_[i]->writeToXmlForReloadParticle(filefullpath);
		}
	}
	//=============================================================================================//
	void ReloadParticleIO::readFromFile(size_t restart_step)
	{
		std::cout << "\n Reloading particles from files." << std::endl;
		for (size_t i = 0; i < bodies_.size(); ++i)
		{
			std::string filefullpath = file_paths_[i];

			if (!fs::exists(filefullpath))
			{
				std::cout << "\n Error: the input file:" << filefullpath << " is not exists" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}

			bodies_[i]->readFromXmlForReloadParticle(filefullpath);
		}
	}
	//=============================================================================================//
	RestartIO::RestartIO(In_Output &in_output, SPHBodyVector bodies)
		: BodyStatesIO(in_output, bodies), overall_file_path_(in_output.restart_folder_ + "/Restart_time_")
	{
		std::transform(bodies.begin(), bodies.end(), std::back_inserter(file_paths_),
					   [&](SPHBody *body) -> std::string
					   { return in_output.restart_folder_ + "/SPHBody_" + body->getBodyName() + "_rst_"; });
	}
	//=============================================================================================//
	void RestartIO::writeToFile(size_t iteration_step)
	{
		std::string overall_filefullpath = overall_file_path_ + std::to_string(iteration_step) + ".dat";
		if (fs::exists(overall_filefullpath))
		{
			fs::remove(overall_filefullpath);
		}
		std::ofstream out_file(overall_filefullpath.c_str(), std::ios::app);
		out_file << std::fixed << std::setprecision(9) << GlobalStaticVariables::physical_time_ << "   \n";
		out_file.close();

		for (size_t i = 0; i < bodies_.size(); ++i)
		{
			std::string filefullpath = file_paths_[i] + std::to_string(iteration_step) + ".xml";

			if (fs::exists(filefullpath))
			{
				fs::remove(filefullpath);
			}
			bodies_[i]->writeParticlesToXmlForRestart(filefullpath);
		}
	}
	//=============================================================================================//
	Real RestartIO::readRestartTime(size_t restart_step)
	{
		std::cout << "\n Reading restart files from the restart step = " << restart_step << std::endl;
		std::string overall_filefullpath = overall_file_path_ + std::to_string(restart_step) + ".dat";
		if (!fs::exists(overall_filefullpath))
		{
			std::cout << "\n Error: the input file:" << overall_filefullpath << " is not exists" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
		Real restart_time;
		std::ifstream in_file(overall_filefullpath.c_str());
		in_file >> restart_time;
		in_file.close();

		return restart_time;
	}
	//=============================================================================================//
	void RestartIO::readFromFile(size_t restart_step)
	{
		for (size_t i = 0; i < bodies_.size(); ++i)
		{
			std::string filefullpath = file_paths_[i] + std::to_string(restart_step) + ".xml";

			if (!fs::exists(filefullpath))
			{
				std::cout << "\n Error: the input file:" << filefullpath << " is not exists" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}

			bodies_[i]->readParticlesFromXmlForRestart(filefullpath);
		}
	}
	//=============================================================================================//
	WriteSimBodyPinData::
		WriteSimBodyPinData(In_Output &in_output, SimTK::Integrator &integ, SimTK::MobilizedBody::Pin &pinbody)
		: WriteSimBodyStates<SimTK::MobilizedBody::Pin>(in_output, integ, pinbody),
		filefullpath_(in_output_.output_folder_ + "/mb_pinbody_data.dat")
	{
		std::ofstream out_file(filefullpath_.c_str(), std::ios::app);

		out_file << "\"time\""
				 << "   ";
		out_file << "  "
				 << "angles"
				 << " ";
		out_file << "  "
				 << "angle_rates"
				 << " ";
		out_file << "\n";

		out_file.close();
	};
	//=============================================================================================//
	void WriteSimBodyPinData::writeToFile(size_t iteration_step)
	{
		std::ofstream out_file(filefullpath_.c_str(), std::ios::app);
		out_file << GlobalStaticVariables::physical_time_ << "   ";
		const SimTK::State &state = integ_.getState();

		out_file << "  " << mobody_.getAngle(state) << "  " << mobody_.getRate(state) << "  ";

		out_file << "\n";
		out_file.close();
	};
	//=============================================================================================//
	WriteSimBodyFreeData::
		WriteSimBodyFreeData(In_Output &in_output,
							 SimTK::Integrator &integ,
							 SimTK::MobilizedBody::Free &pinbody,
							 SimTK::MultibodySystem &mb_system)
		: WriteSimBodyStates<SimTK::MobilizedBody::Free>(in_output, integ, pinbody),
		filefullpath_(in_output_.output_folder_ + "/mb_Freebody_data.dat"),
		mb_system_(mb_system)
	{
		std::ofstream out_file(filefullpath_.c_str(), std::ios::app);
		// header
		out_file << "#THE UNITS PRESENTED BELOW ARE JUST FOR REFERENCE. THE TRUE UNITS ARE " 
				 << "DEPENDENT ON THE UNITS USED TO SET UP THE SIMULATION\n";
		
		out_file << "time [s], "
				 << "COM.x [m], "
				 << "COM.y [m], "
				 << "COM.z [m], "
				 << "vel.x [m/s], "
				 << "vel.y [m/s], "
				 << "vel.z [m/s], "
				 << "acc.x [m/s^2], "
				 << "acc.y [m/s^2], "
				 << "acc.z [m/s^2], "
				 << "rot.roll [rad], "
				 << "rot.pitch [rad], "
				 << "rot.yaw [rad], "
				 << "angvel.x [rad/s], "
				 << "angvel.y [rad/s], "
				 << "angvel.z [rad/s], "
				 << "angacc.x [rad/s^2], "
				 << "angacc.y [rad/s^2], "
				 << "angacc.z [rad/s^2], "
				 << "\n";

		out_file.close();
	};
	//=============================================================================================//
	void WriteSimBodyFreeData::writeToFile(const size_t iteration_step)
	{
		std::ofstream out_file(filefullpath_.c_str(), std::ios::app);

		const SimTK::State &state = integ_.getState();
		mb_system_.realize(state, SimTK::Stage::Acceleration);

		const SimTK::Vec3 &pos = mobody_.getBodyOriginLocation(state);
		const SimTK::SpatialVec &vel = mobody_.getBodyVelocity(state);		// < {angularVel, linearVel}
		const SimTK::SpatialVec &acc = mobody_.getBodyAcceleration(state);  // < {angularAcc, linearAcc}

		const SimTK::Rotation &rot = mobody_.getBodyRotation(state);
		const SimTK::Vec3 angles = rot				// < Vec3(roll, pitch, yaw) 
			.convertThreeAxesRotationToThreeAngles(
				SimTK::BodyOrSpaceType::BodyRotationSequence,
				SimTK::XAxis,
				SimTK::YAxis,
				SimTK::ZAxis
			);

		out_file << GlobalStaticVariables::physical_time_ << ", "
		
				 << pos[0] << ", "
				 << pos[1] << ", "
				 << pos[2] << ", "

				 << vel[1][0] << ", "
				 << vel[1][1] << ", "
				 << vel[1][2] << ", "

				 << acc[1][0] << ", "
				 << acc[1][1] << ", "
				 << acc[1][2] << ", "

				 << angles[0] << ", "
				 << angles[1] << ", "
				 << angles[2] << ", "

				 << vel[0][0] << ", "
				 << vel[0][1] << ", "
				 << vel[0][2] << ", "

				 << acc[0][0] << ", "
				 << acc[0][1] << ", "
				 << acc[0][2] << "\n";

		out_file.close();
	};
	//=================================================================================================//
	ReloadMaterialParameterIO::ReloadMaterialParameterIO(In_Output &in_output, SharedPtr<BaseMaterial> material)
		: in_output_(in_output), material_(material.get()),
		  file_path_(in_output.reload_folder_ + "/Material_" + material->LocalParametersName() + "_rld.xml") {}
	//=================================================================================================//
	ReloadMaterialParameterIO::
		ReloadMaterialParameterIO(In_Output &in_output, SharedPtr<BaseMaterial> material, const std::string &given_parameters_name)
		: in_output_(in_output), material_(material.get()),
		  file_path_(in_output.reload_folder_ + "/Material_" + given_parameters_name + "_rld.xml") {}
	//=================================================================================================//
	void ReloadMaterialParameterIO::writeToFile(size_t iteration_step)
	{
		std::string reload_material_folder = in_output_.reload_folder_;
		if (!fs::exists(reload_material_folder))
		{
			fs::create_directory(reload_material_folder);
		}

		if (fs::exists(file_path_))
		{
			fs::remove(file_path_);
		}
		material_->writeToXmlForReloadLocalParameters(file_path_);
	}
	//=================================================================================================//
	void ReloadMaterialParameterIO::readFromFile(size_t restart_step)
	{
		if (!fs::exists(file_path_))
		{
			std::cout << "\n Error: the reloading material property file:" << file_path_ << " is not exists" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
		material_->readFromXmlForLocalParameters(file_path_);
	}
	//=================================================================================================//
}
