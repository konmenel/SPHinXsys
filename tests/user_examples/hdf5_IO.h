#pragma once
#include <string>
#include <vector>
#include <map>

#include "H5Cpp.h"
#include "H5File.hpp"


/*
NOTE:
Method:
```
std::vector<double> chrono::collision::ChCollisionModelBullet::GetShapeDimensions(int index) const;
```
Returns:
```
SPHERE       radius
BOX          x-halfdim y-halfdim z-halfdim
ELLIPSOID    x-radius y-radius z-radius
CYLINDER     x-radius z-radius halflength
CYLSHELL     radius halflength
```

Collision Shape Types enum:
```
enum chrono::collision::ChCollisionShape::Type {
    SPHERE,
    ELLIPSOID,
    BOX,
    CYLINDER,
    CYLSHELL,
    CONVEXHULL,
    TRIANGLEMESH,
    BARREL,
    POINT,
    TRIANGLE,
    CAPSULE,
    CONE,         // Currently implemented in parallel only
    ROUNDEDBOX,   // Currently implemented in parallel only
    ROUNDEDCYL,   // Currently implemented in parallel only
    ROUNDEDCONE,  // Currently implemented in parallel only
    CONVEX,       // Currently implemented in parallel only
    TETRAHEDRON,  // Currently implemented in parallel only
    PATH2D,
    UNKNOWN_SHAPE
};
```
*/


namespace HF = HighFive;

namespace SPH {
    

    /**
	 * @class SaveStateToH5
	 * @brief Saves the state of the simulation in a h5 file.
	 */
	class SaveStateToH5
	{
	public:
		SaveStateToH5(In_Output &in_output)
        : in_output_(in_output), save_counter_(0),
          bodies_(in_output_.sph_system_.bodies_) 
        {
            stationary_bodies_ = std::vector<bool>(bodies_.size(), false);
        }
		
        ~SaveStateToH5(){};

        void setStationaryBody(size_t index) { stationary_bodies_[index] = true; }
        
        // Note: If two bodies have the same name they will both be marked as stationary
        void setStationaryBody(const std::string &bodyname) 
        {
            for (size_t i = 0; i < bodies_.size(); i++) {
                if (bodyname == bodies_[i]->getBodyName()) {
                    stationary_bodies_[i] = true; 
                }
            }
        }

        void setMovingBody(size_t index) { stationary_bodies_[index] = false; }
        
        // Note: If two bodies have the same name they will both be marked as moving
        void setMovingBody(const std::string &bodyname) 
        {
            for (size_t i = 0; i < bodies_.size(); i++) {
                if (bodyname == bodies_[i]->getBodyName()) {
                    stationary_bodies_[i] = false; 
                }
            }
        }

        void writeConfiguration()
        {
            using std::string;
            using std::vector;
            using namespace chrono::collision;
            
            string filefullpath = in_output_.output_folder_ + "/Config.h5";
            HF::File out_file(filefullpath, HF::File::ReadWrite | HF::File::Create | HF::File::Truncate);

            // Either 2 or 3 for 2D or 3D
            int sim_type = Vecd::size();
            HF::DataSet sim_type_attr = out_file.createDataSet <int>("SimulationType", HF::DataSpace::From(sim_type));
            sim_type_attr.write(sim_type);

            // Precision used in simulation (1 if Real=float or 2 if Real=double)
            int precision = SimTK_DEFAULT_PRECISION;
            HF::DataSet precision_attr = out_file.createDataSet<int>("Precision", HF::DataSpace::From(precision));
            precision_attr.write(precision);

            // Chrono shape data for chrono bodies
            HF::Group chrono_bodies_group = out_file.createGroup("ChronoShapes");
            for (size_t i = 0; i < bodies_.size(); i++) {
                SolidBodyWithChrono *body = dynamic_cast<SolidBodyWithChrono*>(bodies_[i]);
                
                if (body) {
                    HF::Group body_group = chrono_bodies_group.createGroup(std::to_string(i));
                    auto collision_model = body->getChBody()->GetCollisionModel();  // < std::shared_ptr
                    const auto &shapes = collision_model->GetShapes();              // < std::vector
                    for (int j = 0; j < shapes.size(); j++) {
                        int type = shapes[j]->GetType();
                        HF::DataSet type_dataset = body_group.createDataSet<int>("ShapeType", HF::DataSpace::From(type));
                        type_dataset.write(type);

                        std::vector<Real> shape_dims = collision_model->GetShapeDimensions(j);
                        HF::DataSpace data_dims = { shape_dims.size() };
                        HF::DataSet shape_dims_dataset = body_group.createDataSet<int>("ShapeDims", data_dims);
                        shape_dims_dataset.write(shape_dims);
                    }
                }
            }
        }

        void writeToFile()
        {
            using std::string;
            using std::vector;

            string filefullpath = in_output_.output_folder_ + "/State_" + std::to_string(save_counter_) + ".h5";
            HF::File out_file(filefullpath, HF::File::ReadWrite | HF::File::Create | HF::File::Truncate);

            // Write time
            HF::DataSet time_dataset = out_file.createDataSet<Real>("Time",
                HF::DataSpace::From(GlobalStaticVariables::physical_time_));
            time_dataset.write(GlobalStaticVariables::physical_time_);

            // Creating the three groups.
            HF::Group fluids_group = out_file.createGroup("Fluids");
            HF::Group solids_group = out_file.createGroup("Solids");
            HF::Group fictitious_group = out_file.createGroup("Fictitious");

            for (size_t i = 0; i < bodies_.size(); i++) {
                if (!stationary_bodies_[i]) {    // Case of moving body
                    writeMovingBody(out_file, i);
                
                } else {    // Else stationary case
                    if (save_counter_ != 0) {
                        writeStationaryBody(out_file, i);
                    } else {
                        writeMovingBody(out_file, i);
                    }
                }
            }
            save_counter_++;
        }


	private:
        In_Output &in_output_;
        size_t save_counter_;
        SPHBodyVector &bodies_;
        std::vector<bool> stationary_bodies_;

        void writeMovingBody(HF::File &out_file, size_t index)
        {
            SPHBody *body = bodies_[index];
            BaseParticles *particles = body->base_particles_;

            HF::Group body_group;
            string body_name = body->getBodyName();
            
            // Check body type
            if (dynamic_cast<FluidBody *>(body)) {
                body_group = out_file.getGroup("Fluids").createGroup(std::to_string(index));
            
            } else if (dynamic_cast<SolidBody *>(body)) {
                body_group = out_file.getGroup("Solids").createGroup(std::to_string(index));

                // Write Chrono data
                if (auto chrono_body = dynamic_cast<SolidBodyWithChrono *>(body)) {
                    HF::Group chrono_group = body_group.createGroup("ChronoData");
                    chrono_body->writeChronoDataToH5(chrono_group);
                }

            } else if (dynamic_cast<FictitiousBody *>(body)) {
                body_group = out_file.getGroup("Fictitious").createGroup(std::to_string(index));
            
            } else {
                std::cout << "[WARNING] Unhandled body during save state with name: " << body_name
                << "\nNo match found for the body. Body data not saved!" << std::endl;
            }
            
            HF::Attribute name_attr = body_group.createAttribute<string>("Name", HF::DataSpace::From(body_name));
            name_attr.write(body_name);

            particles->writeStateToH5(body_group);
        }

        void writeStationaryBody(HF::File &out_file, size_t index)
        {
            SPHBody *body = bodies_[index];
            std::string group_name;
            std::string index_str = std::to_string(index);

            if (dynamic_cast<FluidBody *>(body)) {
                group_name = "Fluid";
                
            } else if (dynamic_cast<SolidBody *>(body)) {
                group_name = "Solids";


            } else if (dynamic_cast<FictitiousBody *>(body)) {
                group_name = "Fictitious";
            
            } else {
                std::cout << "[WARNING] Unhandled body during save state with name: " << body->getBodyName()
                << "\nNo match found for the body. Body data not saved!" << std::endl;
            }

            
            out_file.getGroup(group_name).createExternalLink(index_str, "State_0.h5", 
                "/" + group_name + "/" + index_str);
        }
	};
}
