/**
 * @file base_particles.cpp
 * @brief Definition of functions declared in base_particles.h
 * @author	Xiangyu Hu and Chi Zhang
 */
#include "base_particles.h"
#include "base_particles.hpp"

#include "base_body.h"
#include "base_material.h"
#include "base_particle_generator.h"
#include "xml_engine.h"

namespace SPH
{
	//=================================================================================================//
	BaseParticles::BaseParticles(SPHBody &sph_body,
								 SharedPtr<BaseMaterial> base_material_ptr,
								 SharedPtr<ParticleGenerator> particle_generator)
		: base_material_(base_material_ptr_keeper_.assignPtr(base_material_ptr)),
		  rho0_(base_material_->ReferenceDensity()),
		  sigma0_(sph_body.sph_adaptation_->ReferenceNumberDensity()),
		  speed_max_(0.0), signal_speed_max_(0.0),
		  total_real_particles_(0), real_particles_bound_(0), total_ghost_particles_(0),
		  sph_body_(&sph_body), body_name_(sph_body.getBodyName()),
		  restart_xml_engine_("xml_restart", "particles"),
		  reload_xml_engine_("xml_particle_reload", "particles")
	{
		sph_body.assignBaseParticles(this);
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable<indexVector, Vecd>(pos_n_, "Position");
		registerAVariable<indexVector, Vecd>(vel_n_, "Velocity");
		registerAVariable<indexVector, Vecd>(dvel_dt_, "Acceleration");
		registerAVariable<indexVector, Vecd>(dvel_dt_prior_, "PriorAcceleration");
		registerAVariable<indexScalar, Real>(Vol_, "Volume");
		registerAVariable<indexScalar, Real>(rho_n_, "Density");
		registerAVariable<indexScalar, Real>(mass_, "Mass");
		//----------------------------------------------------------------------
		//		add basic output particle data
		//----------------------------------------------------------------------
		addAVariableToWrite<indexVector, Vecd>("Velocity");
		addAVariableToWrite<indexVector, Vecd>("Acceleration");
		//----------------------------------------------------------------------
		//		add restart output particle data
		//----------------------------------------------------------------------
		addAVariableNameToList<indexVector, Vecd>(variables_to_restart_, "Position");
		addAVariableNameToList<indexVector, Vecd>(variables_to_restart_, "Velocity");
		addAVariableNameToList<indexVector, Vecd>(variables_to_restart_, "Acceleration");
		addAVariableNameToList<indexScalar, Real>(variables_to_restart_, "Volume");

		particle_generator->initialize(&sph_body);
		particle_generator->createBaseParticles(this);
		real_particles_bound_ = total_real_particles_;

		sph_body.sph_adaptation_->assignBaseParticles(this);
		base_material_->assignBaseParticles(this);
	}
	//=================================================================================================//
	void BaseParticles::initializeABaseParticle(Vecd pnt, Real Vol_0)
	{
		total_real_particles_++;
		sequence_.push_back(0);
		sorted_id_.push_back(pos_n_.size());
		unsorted_id_.push_back(pos_n_.size());

		pos_n_.push_back(pnt);
		vel_n_.push_back(Vecd(0));
		dvel_dt_.push_back(Vecd(0));
		dvel_dt_prior_.push_back(Vecd(0));

		Vol_.push_back(Vol_0);
		rho_n_.push_back(rho0_);
		mass_.push_back(rho0_ * Vol_0);
	}
	//=================================================================================================//
	void BaseParticles::addAParticleEntry()
	{
		sequence_.push_back(0);
		sorted_id_.push_back(pos_n_.size());
		unsorted_id_.push_back(pos_n_.size());

		add_a_particle_value_(all_particle_data_);
	}
	//=================================================================================================//
	void BaseParticles::addBufferParticles(size_t buffer_size)
	{
		for (size_t i = 0; i != buffer_size; ++i)
		{
			addAParticleEntry();
		}
		real_particles_bound_ += buffer_size;
	}
	//=================================================================================================//
	void BaseParticles::copyFromAnotherParticle(size_t this_index, size_t another_index)
	{
		updateFromAnotherParticle(this_index, another_index);
	}
	//=================================================================================================//
	void BaseParticles::updateFromAnotherParticle(size_t this_index, size_t another_index)
	{
		copy_a_particle_value_(all_particle_data_, this_index, another_index);
	}
	//=================================================================================================//
	size_t BaseParticles::insertAGhostParticle(size_t index_i)
	{
		total_ghost_particles_ += 1;
		size_t expected_size = real_particles_bound_ + total_ghost_particles_;
		size_t expected_particle_index = expected_size - 1;
		if (expected_size <= pos_n_.size())
		{
			copyFromAnotherParticle(expected_particle_index, index_i);
			/** For a ghost particle, its sorted id is that of corresponding real particle. */
			sorted_id_[expected_particle_index] = index_i;
		}
		else
		{
			addAParticleEntry();
			copyFromAnotherParticle(expected_particle_index, index_i);
			/** For a ghost particle, its sorted id is that of corresponding real particle. */
			sorted_id_[expected_particle_index] = index_i;
		}
		return expected_particle_index;
	}
	//=================================================================================================//
	void BaseParticles::switchToBufferParticle(size_t index_i)
	{
		size_t last_real_particle_index = total_real_particles_ - 1;
		updateFromAnotherParticle(index_i, last_real_particle_index);
		unsorted_id_[index_i] = unsorted_id_[last_real_particle_index];
		total_real_particles_ -= 1;
	}
//=================================================================================================//
	void BaseParticles::writeParticlesToVtuFile(std::ostream& output_file)
	{
		size_t total_real_particles = total_real_particles_;

		//write current/final particle positions first
		output_file << "   <Points>\n";
		output_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_real_particles; ++i) {
			Vec3d particle_position = upgradeToVector3D(pos_n_[i]);
			output_file << particle_position[0] << " " << particle_position[1] << " " << particle_position[2] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
		output_file << "   </Points>\n";

		//write header of particles data
		output_file << "   <PointData  Vectors=\"vector\">\n";

		//write sorted particles ID
		output_file << "    <DataArray Name=\"SortedParticle_ID\" type=\"Int32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_real_particles; ++i) {
			output_file << i << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		//write unsorted particles ID
		output_file << "    <DataArray Name=\"UnsortedParticle_ID\" type=\"Int32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_real_particles; ++i) {
			output_file << unsorted_id_[i] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		//write matrices
		for (std::pair<std::string, size_t>& name_index : variables_to_write_[indexMatrix])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Matd>& variable = *(std::get<indexMatrix>(all_particle_data_)[name_index.second]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Float32\"  NumberOfComponents=\"9\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_real_particles; ++i) {
				Mat3d matrix_value = upgradeToMatrix3D(variable[i]);
				for (int k = 0; k != 3; ++k) {
					Vec3d col_vector = matrix_value.col(k);
					output_file << std::fixed << std::setprecision(9) << col_vector[0] << " " << col_vector[1] << " " << col_vector[2] << " ";
				}
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}

		//write vectors
		for (std::pair<std::string, size_t>& name_index : variables_to_write_[indexVector])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Vecd>& variable = *(std::get<indexVector>(all_particle_data_)[name_index.second]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_real_particles; ++i) {
				Vec3d vector_value = upgradeToVector3D(variable[i]);
				output_file << std::fixed << std::setprecision(9) << vector_value[0] << " " << vector_value[1] << " " << vector_value[2] << " ";
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}

		//write scalars
		for (std::pair<std::string, size_t>& name_index : variables_to_write_[indexScalar])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Real>& variable = *(std::get<indexScalar>(all_particle_data_)[name_index.second]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Float32\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_real_particles; ++i) {
				output_file << std::fixed << std::setprecision(9) << variable[i] << " ";
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}

		//write integers
		for (std::pair<std::string, size_t>& name_index : variables_to_write_[indexInteger])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<int>& variable = *(std::get<indexInteger>(all_particle_data_)[name_index.second]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Int32\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_real_particles; ++i) {
				output_file << std::fixed << std::setprecision(9) << variable[i] << " ";
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}
	}
	//=================================================================================================//
	void BaseParticles::writeParticlesToVtpFile(std::ofstream &output_file)
	{
		size_t total_real_particles = total_real_particles_;

		//write current/final particle positions first
		output_file << "   <Points>\n";
		output_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_real_particles; ++i)
		{
			Vec3d particle_position = upgradeToVector3D(pos_n_[i]);
			output_file << particle_position[0] << " " << particle_position[1] << " " << particle_position[2] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
		output_file << "   </Points>\n";

		//write header of particles data
		output_file << "   <PointData  Vectors=\"vector\">\n";

		//write sorted particles ID
		output_file << "    <DataArray Name=\"SortedParticle_ID\" type=\"Int32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_real_particles; ++i)
		{
			output_file << i << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		//write unsorted particles ID
		output_file << "    <DataArray Name=\"UnsortedParticle_ID\" type=\"Int32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_real_particles; ++i)
		{
			output_file << unsorted_id_[i] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		//write matrices
		for (std::pair<std::string, size_t> &name_index : variables_to_write_[indexMatrix])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Matd> &variable = *(std::get<indexMatrix>(all_particle_data_)[name_index.second]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Float32\"  NumberOfComponents=\"9\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_real_particles; ++i)
			{
				Mat3d matrix_value = upgradeToMatrix3D(variable[i]);
				for (int k = 0; k != 3; ++k)
				{
					Vec3d col_vector = matrix_value.col(k);
					output_file << std::fixed << std::setprecision(9) << col_vector[0] << " " << col_vector[1] << " " << col_vector[2] << " ";
				}
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}

		//write vectors
		for (std::pair<std::string, size_t> &name_index : variables_to_write_[indexVector])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Vecd> &variable = *(std::get<indexVector>(all_particle_data_)[name_index.second]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_real_particles; ++i)
			{
				Vec3d vector_value = upgradeToVector3D(variable[i]);
				output_file << std::fixed << std::setprecision(9) << vector_value[0] << " " << vector_value[1] << " " << vector_value[2] << " ";
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}

		//write scalars
		for (std::pair<std::string, size_t> &name_index : variables_to_write_[indexScalar])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Real> &variable = *(std::get<indexScalar>(all_particle_data_)[name_index.second]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Float32\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_real_particles; ++i)
			{
				output_file << std::fixed << std::setprecision(9) << variable[i] << " ";
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}

		//write integers
		for (std::pair<std::string, size_t> &name_index : variables_to_write_[indexInteger])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<int> &variable = *(std::get<indexInteger>(all_particle_data_)[name_index.second]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Int32\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_real_particles; ++i)
			{
				output_file << std::fixed << std::setprecision(9) << variable[i] << " ";
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}
	}
	//=================================================================================================//
	void BaseParticles::writeParticlesToBinLecVtkFile(std::ofstream &output_file)
	{
		// Getting the endianness of the machine
		const static Endianness endianness = Endian::getSystemEndianness();

		// Create the buffers
		const int BUFF_CAP = 1024;
		float bufferf[BUFF_CAP];
		int32_t bufferi[BUFF_CAP];
		uint32_t bufferui[BUFF_CAP];
		size_t bufferf_size = 0;
		size_t bufferi_size = 0;
		size_t bufferui_size = 0;
		
		size_t total_real_particles = total_real_particles_;

		//function to flush the buffer
		auto flush_buffer = [&output_file](void *buf, size_t bytes_count, size_t type_size)
		{
			if (endianness == Endianness::little) {
				Endian::writeDataReverseEndianness(output_file, buf, bytes_count, type_size);
			} else {
				output_file.write(reinterpret_cast<const char *>(buf), bytes_count);
			}
		};

		//write current/final particle positions first
		output_file << "POINTS " << total_real_particles << " float\n";
		for (size_t i = 0; i != total_real_particles; ++i)
		{
			Vec3d particle_position = upgradeToVector3D(pos_n_[i]);
			for (size_t j = 0; j < 3; ++j)
			{
				bufferf[bufferf_size++] = static_cast<float>(particle_position[j]);

				if (bufferf_size == BUFF_CAP) {
					flush_buffer(bufferf, bufferf_size * sizeof(float), sizeof(float));
					bufferf_size = 0;
				}
				// if (endianness == Endianness::little) {
				// 	Endian::writeDataReverseEndianness(output_file, &part_pos_float, sizeof(float), sizeof(float));
				// } else {
				// 	output_file.write(reinterpret_cast<const char *>(&part_pos_float), sizeof(float));
				// }
			}
		}
		if (bufferf_size > 0) {
			flush_buffer(bufferf, bufferf_size * sizeof(float), sizeof(float));
			bufferf_size = 0;
		}
		output_file << "\n";

		//write the vertices
		const int32_t number_of_points = 1;
		output_file << "VERTICES " << total_real_particles << " " << total_real_particles * 2 << "\n";
		for (int32_t i = 0; i != total_real_particles; ++i)
		{
			bufferi[bufferi_size++] = number_of_points;

			if  (bufferi_size == BUFF_CAP) {
				flush_buffer(bufferi, bufferi_size * sizeof(int32_t), sizeof(int32_t));
				bufferi_size = 0;
			}

			bufferi[bufferi_size++] = i;
			if  (bufferi_size == BUFF_CAP) {
				flush_buffer(bufferi, bufferi_size * sizeof(int32_t), sizeof(int32_t));
				bufferi_size = 0;
			}
			// Endian::writeDataReverseEndianness(output_file, &number_of_points, sizeof(int), sizeof(int));
			// Endian::writeDataReverseEndianness(output_file, &i, sizeof(int), sizeof(int));
		}
		if (bufferi_size > 0) {
			flush_buffer(bufferi, bufferi_size * sizeof(int32_t), sizeof(int32_t));
			bufferi_size = 0;
		}
		output_file << "\n";


		//write sorted particles ID
		output_file << "POINT_DATA " << total_real_particles << "\n";
		output_file << "SCALARS Particle_ID unsigned_int\n";
		output_file << "LOOKUP_TABLE default\n";
		for (uint32_t i = 0; i != total_real_particles; ++i)
		{
			bufferui[bufferui_size++] = i;

			if  (bufferui_size == BUFF_CAP) {
				flush_buffer(bufferui, bufferui_size * sizeof(uint32_t), sizeof(uint32_t));
				bufferui_size = 0;
			}
			// if (endianness == Endianness::little) {
			// 	Endian::writeDataReverseEndianness(output_file, &i,
			// 		sizeof(uint32_t), sizeof(uint32_t));
			// } else {
			// 	output_file.write(reinterpret_cast<const char *>(&i), sizeof(uint32_t));
			// }
		}
		if (bufferui_size > 0) {
			flush_buffer(bufferui, bufferui_size * sizeof(uint32_t), sizeof(uint32_t));
			bufferui_size = 0;
		}
		output_file << "\n";

		//write unsorted particles ID
		output_file << "SCALARS UnsortedParticle_ID unsigned_int\n";
		output_file << "LOOKUP_TABLE default\n";
		for (size_t i = 0; i != total_real_particles; ++i)
		{
			bufferui[bufferui_size++] = unsorted_id_[i];

			if  (bufferui_size == BUFF_CAP) {
				flush_buffer(bufferui, bufferui_size * sizeof(uint32_t), sizeof(uint32_t));
				bufferui_size = 0;
			}
			// uint32_t id = unsorted_id_[i];
			// if (endianness == Endianness::little) {
			// 	Endian::writeDataReverseEndianness(output_file, 
			// 		reinterpret_cast<const char *>(&id),
			// 		sizeof(uint32_t), sizeof(uint32_t));
			// } else {
			// 	output_file.write(reinterpret_cast<const char *>(&id), sizeof(uint32_t));
			// }
		}
		if (bufferui_size > 0) {
			flush_buffer(bufferui, bufferui_size * sizeof(uint32_t), sizeof(uint32_t));
			bufferui_size = 0;
		}
		output_file << "\n";

		//write header of field data
		int fields_size = variables_to_write_[indexMatrix].size()
			+ variables_to_write_[indexVector].size()
			+ variables_to_write_[indexScalar].size()
			+ variables_to_write_[indexInteger].size();
		output_file << "FIELD FieldData " << fields_size << "\n";

		//write matrices
		for (std::pair<std::string, size_t> &name_index : variables_to_write_[indexMatrix])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Matd> &variable = *(std::get<indexMatrix>(all_particle_data_)[name_index.second]);

			// write header
			output_file << variable_name << " 9 " << total_real_particles << " float\n";
			
			for (size_t i = 0; i != total_real_particles; ++i)
			{
				Mat3d matrix_value = upgradeToMatrix3D(variable[i]);
				for (int k = 0; k != 3; ++k)
				{
					Vec3d col_vector = matrix_value.col(k);

					for (int j = 0; j < 3; ++j)
					{
						bufferf[bufferf_size++] = static_cast<float>(col_vector[j]);

						if  (bufferf_size == BUFF_CAP) {
							flush_buffer(bufferf, bufferf_size * sizeof(float), sizeof(float));
							bufferf_size = 0;
						}

						// if (endianness == Endianness::little) {
						// 	Endian::writeDataReverseEndianness(output_file, &elem_float,
						// 		sizeof(float), sizeof(float));
						// } else {
						// 	output_file.write(reinterpret_cast<const char *>(&elem_float), sizeof(float));
						// }
					}
				}
			}
			if (bufferf_size > 0) {
				flush_buffer(bufferf, bufferf_size * sizeof(float), sizeof(float));
				bufferf_size = 0;
			}
			output_file << "\n";
		}

		//write vectors
		for (std::pair<std::string, size_t> &name_index : variables_to_write_[indexVector])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Vecd> &variable = *(std::get<indexVector>(all_particle_data_)[name_index.second]);
			
			// write header
			output_file << variable_name << " 3 " << total_real_particles << " float\n";

			for (size_t i = 0; i != total_real_particles; ++i)
			{
				Vec3d vector_value = upgradeToVector3D(variable[i]);
				
				for (int j = 0; j < 3; ++j)
				{
					bufferf[bufferf_size++] = static_cast<float>(vector_value[j]);

					if  (bufferf_size == BUFF_CAP) {
						flush_buffer(bufferf, bufferf_size * sizeof(float), sizeof(float));
						bufferf_size = 0;
					}

					// if (endianness == Endianness::little) {
					// 	Endian::writeDataReverseEndianness(output_file, &elem_float,
					// 		sizeof(float), sizeof(float));
					// } else {
					// 	output_file.write(reinterpret_cast<const char *>(&elem_float), sizeof(float));
					// }
				}
			}
			if (bufferf_size > 0) {
				flush_buffer(bufferf, bufferf_size * sizeof(float), sizeof(float));
				bufferf_size = 0;
			}
			output_file << "\n";
		}

		//write scalars
		for (std::pair<std::string, size_t> &name_index : variables_to_write_[indexScalar])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Real> &variable = *(std::get<indexScalar>(all_particle_data_)[name_index.second]);

			// write header
			output_file << variable_name << " 1 " << total_real_particles << " float\n";

			for (size_t i = 0; i != total_real_particles; ++i)
			{
				bufferf[bufferf_size++] = static_cast<float>(variable[i]);

				if  (bufferf_size == BUFF_CAP) {
					flush_buffer(bufferf, bufferf_size * sizeof(float), sizeof(float));
					bufferf_size = 0;
				}

				// if (endianness == Endianness::little) {
				// 	Endian::writeDataReverseEndianness(output_file, &scalar_value,
				// 		sizeof(float), sizeof(float));
				// } else {
				// 	output_file.write(reinterpret_cast<const char *>(&scalar_value), sizeof(float));
				// }
			}
			if (bufferf_size > 0) {
				flush_buffer(bufferf, bufferf_size * sizeof(float), sizeof(float));
				bufferf_size = 0;
			}
			output_file << "\n";
		}

		//write integers
		for (std::pair<std::string, size_t> &name_index : variables_to_write_[indexInteger])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<int> &variable = *(std::get<indexInteger>(all_particle_data_)[name_index.second]);
			
			// write header
			output_file << variable_name << " 1 " << total_real_particles << " int\n";

			for (size_t i = 0; i != total_real_particles; ++i)
			{
				// make sure correct number of bytes are written
				bufferi[bufferi_size++] = variable[i];

				if  (bufferi_size == BUFF_CAP) {
					flush_buffer(bufferi, bufferi_size * sizeof(int), sizeof(int));
					bufferi_size = 0;
				}

				// if (endianness == Endianness::little) {
				// 	Endian::writeDataReverseEndianness(output_file, &integer_value,
				// 		sizeof(int), sizeof(int));
				// } else {
				// 	output_file.write(reinterpret_cast<const char *>(&integer_value), sizeof(int));
				// }
			}
			if (bufferi_size > 0) {
				flush_buffer(bufferi, bufferi_size * sizeof(int), sizeof(int));
				bufferi_size = 0;
			}
			output_file << "\n";
		}
	}
	//=================================================================================================//
	void BaseParticles::writePltFileHeader(std::ofstream &output_file)
	{
		output_file << " VARIABLES = \"x\",\"y\",\"z\",\"ID\"";

		for (size_t l = 0; l != variables_to_write_[indexInteger].size(); ++l)
		{
			std::string variable_name = variables_to_write_[indexInteger][l].first;
			output_file << ",\"" << variable_name << "\"";
		};

		for (size_t l = 0; l != variables_to_write_[indexVector].size(); ++l)
		{
			std::string variable_name = variables_to_write_[indexVector][l].first;
			output_file << ",\"" << variable_name << "_x\""
						<< ",\"" << variable_name << "_y\""
						<< ",\"" << variable_name << "_z\"";
		};
		for (size_t l = 0; l != variables_to_write_[indexScalar].size(); ++l)
		{
			std::string variable_name = variables_to_write_[indexScalar][l].first;
			output_file << ",\"" << variable_name << "\"";
		};
	}
	//=================================================================================================//
	void BaseParticles::writePltFileParticleData(std::ofstream &output_file, size_t index_i)
	{
		//write particle positions and index first
		Vec3d particle_position = upgradeToVector3D(pos_n_[index_i]);
		output_file << particle_position[0] << " " << particle_position[1] << " " << particle_position[2] << " "
					<< index_i << " ";

		for (std::pair<std::string, size_t> &name_index : variables_to_write_[indexInteger])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<int> &variable = *(std::get<indexInteger>(all_particle_data_)[name_index.second]);
			output_file << variable[index_i] << " ";
		};

		for (std::pair<std::string, size_t> &name_index : variables_to_write_[indexVector])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Vecd> &variable = *(std::get<indexVector>(all_particle_data_)[name_index.second]);
			Vec3d vector_value = upgradeToVector3D(variable[index_i]);
			output_file << vector_value[0] << " " << vector_value[1] << " " << vector_value[2] << " ";
		};

		for (std::pair<std::string, size_t> &name_index : variables_to_write_[indexScalar])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Real> &variable = *(std::get<indexScalar>(all_particle_data_)[name_index.second]);
			output_file << variable[index_i] << " ";
		};
	}
	//=================================================================================================//
	void BaseParticles::writeParticlesToPltFile(std::ofstream &output_file)
	{
		writePltFileHeader(output_file);
		output_file << "\n";

		size_t total_real_particles = total_real_particles_;
		for (size_t i = 0; i != total_real_particles; ++i)
		{
			writePltFileParticleData(output_file, i);
			output_file << "\n";
		};
	}
	//=================================================================================================//
	void BaseParticles::writeSurfaceParticlesToVtuFile(std::ofstream& output_file, BodySurface& surface_particles)
	{
		size_t total_surface_particles = surface_particles.body_part_particles_.size();

		//write current/final particle positions first
		output_file << "   <Points>\n";
		output_file << "    <DataArray Name=\"Position\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_surface_particles; ++i) {
			size_t particle_i = surface_particles.body_part_particles_[i];
			Vec3d particle_position = upgradeToVector3D(pos_n_[particle_i]);
			output_file << particle_position[0] << " " << particle_position[1] << " " << particle_position[2] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
		output_file << "   </Points>\n";

		//write header of particles data
		output_file << "   <PointData  Vectors=\"vector\">\n";

		//write sorted particles ID
		output_file << "    <DataArray Name=\"SortedParticle_ID\" type=\"Int32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_surface_particles; ++i) {
			size_t particle_i = surface_particles.body_part_particles_[i];
			output_file << particle_i << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		//write unsorted particles ID
		output_file << "    <DataArray Name=\"UnsortedParticle_ID\" type=\"Int32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_surface_particles; ++i) {
			size_t particle_i = surface_particles.body_part_particles_[i];
			output_file << unsorted_id_[particle_i] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		//write matrices
		for (std::pair<std::string, size_t>& name_index : variables_to_write_[indexMatrix])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Matd>& variable = *(std::get<indexMatrix>(all_particle_data_)[name_index.second]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Float32\"  NumberOfComponents=\"9\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_surface_particles; ++i) {
				size_t particle_i = surface_particles.body_part_particles_[i];
				Mat3d matrix_value = upgradeToMatrix3D(variable[particle_i]);
				for (int k = 0; k != 3; ++k) {
					Vec3d col_vector = matrix_value.col(k);
					output_file << std::fixed << std::setprecision(9) << col_vector[0] << " " << col_vector[1] << " " << col_vector[2] << " ";
				}
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}

		//write vectors
		for (std::pair<std::string, size_t>& name_index : variables_to_write_[indexVector])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Vecd>& variable = *(std::get<indexVector>(all_particle_data_)[name_index.second]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Float32\"  NumberOfComponents=\"3\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_surface_particles; ++i) {
				size_t particle_i = surface_particles.body_part_particles_[i];
				Vec3d vector_value = upgradeToVector3D(variable[particle_i]);
				output_file << std::fixed << std::setprecision(9) << vector_value[0] << " " << vector_value[1] << " " << vector_value[2] << " ";
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}

		//write scalars
		for (std::pair<std::string, size_t>& name_index : variables_to_write_[indexScalar])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<Real>& variable = *(std::get<indexScalar>(all_particle_data_)[name_index.second]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Float32\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_surface_particles; ++i) {
				size_t particle_i = surface_particles.body_part_particles_[i];
				output_file << std::fixed << std::setprecision(9) << variable[particle_i] << " ";
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}

		//write integers
		for (std::pair<std::string, size_t>& name_index : variables_to_write_[indexInteger])
		{
			std::string variable_name = name_index.first;
			StdLargeVec<int>& variable = *(std::get<indexInteger>(all_particle_data_)[name_index.second]);
			output_file << "    <DataArray Name=\"" << variable_name << "\" type=\"Int32\" Format=\"ascii\">\n";
			output_file << "    ";
			for (size_t i = 0; i != total_surface_particles; ++i) {
				size_t particle_i = surface_particles.body_part_particles_[i];
				output_file << std::fixed << std::setprecision(9) << variable[particle_i] << " ";
			}
			output_file << std::endl;
			output_file << "    </DataArray>\n";
		}
	}
	//=================================================================================================//
	void BaseParticles::resizeXmlDocForParticles(XmlEngine &xml_engine)
	{
		size_t total_elements = xml_engine.SizeOfXmlDoc();

		if (total_elements <= total_real_particles_)
		{
			for (size_t i = total_elements; i != total_real_particles_; ++i)
				xml_engine.addElementToXmlDoc("particle");
		}
	}
	//=================================================================================================//
	void BaseParticles::writeParticlesToXmlForRestart(std::string &filefullpath)
	{
		resizeXmlDocForParticles(restart_xml_engine_);
		WriteAParticleVariableToXml write_variable_to_xml(restart_xml_engine_, total_real_particles_);
		ParticleDataOperation<loopVariabaleNameList> loop_variable_namelist;
		loop_variable_namelist(all_particle_data_, variables_to_restart_, write_variable_to_xml);
		restart_xml_engine_.writeToXmlFile(filefullpath);
	}
	//=================================================================================================//
	void BaseParticles::readParticleFromXmlForRestart(std::string &filefullpath)
	{
		restart_xml_engine_.loadXmlFile(filefullpath);
		ReadAParticleVariableFromXml read_variable_from_xml(restart_xml_engine_, total_real_particles_);
		ParticleDataOperation<loopVariabaleNameList> loop_variable_namelist;
		loop_variable_namelist(all_particle_data_, variables_to_restart_, read_variable_from_xml);
	}
	//=================================================================================================//
	void BaseParticles::writeToXmlForReloadParticle(std::string &filefullpath)
	{
		resizeXmlDocForParticles(reload_xml_engine_);
		SimTK::Xml::element_iterator ele_ite = reload_xml_engine_.root_element_.element_begin();
		for (size_t i = 0; i != total_real_particles_; ++i)
		{
			reload_xml_engine_.setAttributeToElement(ele_ite, "Position", pos_n_[i]);
			reload_xml_engine_.setAttributeToElement(ele_ite, "Volume", Vol_[i]);
			ele_ite++;
		}
		reload_xml_engine_.writeToXmlFile(filefullpath);
	}
	//=================================================================================================//
	void BaseParticles::readFromXmlForReloadParticle(std::string &filefullpath)
	{
		reload_xml_engine_.loadXmlFile(filefullpath);
		SimTK::Xml::element_iterator ele_ite = reload_xml_engine_.root_element_.element_begin();
		for (size_t i = 0; i != total_real_particles_; ++i)
		{
			reload_xml_engine_.getRequiredAttributeValue<Vecd>(ele_ite, "Position", pos_n_[i]);
			reload_xml_engine_.getRequiredAttributeValue<Real>(ele_ite, "Volume", Vol_[i]);
			ele_ite++;
		}

		if (reload_xml_engine_.SizeOfXmlDoc() != total_real_particles_)
		{
			std::cout << "\n Error: reload particle number does not match!" << std::endl;
			std::cout << __FILE__ << ':' << __LINE__ << std::endl;
			exit(1);
		}
	}
	//=================================================================================================//
}
