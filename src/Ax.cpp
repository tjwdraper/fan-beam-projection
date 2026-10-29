#include <fstream>
#include <string>
#include <stdexcept>
#include <chrono>
#include <vector>
#include <map>

#include "coord2d.hpp"
#include "json.hpp"
#include "Field.hpp"
#include "fan_beam.hpp"
#include "fan_beam_configuration.hpp"

///////////////////////////////////////////////////////////////////////////////////////////////////
// Read json configuration file and store in structure
///////////////////////////////////////////////////////////////////////////////////////////////////
struct json_config {
    // Path to files
    std::string path_image;
    std::string path_projections;

    // Fan-beam configuration parameters
    fan_beam_configuration geometry;
};

json_config load_config(const std::string& filename) {
    // Open .json configuration file
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Could not open json configuration file: " + filename);
    nlohmann::json json;
    file >> json;

    // Create json configuration structure
    json_config config;

    // Set the image path
    config.path_image = json.at("path_image").get<std::string>();

    // Set the output path
    config.path_projections = json.at("path_projections").get<std::string>();

    // Set fan beam geometry parameters
    if (json.contains("geometry")) {
        const auto& geometry = json.at("geometry");

        if (geometry.contains("DSO"))
            config.geometry.DSO = geometry.at("DSO").get<double>();
        else
            throw std::runtime_error("In json_config load_config(const std::string&) could not find DSO variable in geometry.");

        if (geometry.contains("DSD"))
            config.geometry.DSD = geometry.at("DSD").get<double>();
        else
            throw std::runtime_error("In json_config load_config(const std::string&) could not find DSD variable in geometry.");

        if (geometry.contains("nVoxel")) {
            const std::vector<int> nVoxel = geometry.at("nVoxel").get<std::vector<int>>();
            config.geometry.nVoxel.x = static_cast<unsigned int>(nVoxel[0]);
            config.geometry.nVoxel.y = static_cast<unsigned int>(nVoxel[1]);
        }
        else
            throw std::runtime_error("In json_config load_config(const std::string&) could not find nVoxel variable in geometry.");

        if (geometry.contains("sVoxel")) {
            const std::vector<double> sVoxel = geometry.at("sVoxel").get<std::vector<double>>();
            config.geometry.sVoxel.x = sVoxel[0];
            config.geometry.sVoxel.y = sVoxel[1];
        }
        else
            throw std::runtime_error("In json_config load_config(const std::string&) could not find sVoxel variable in geometry.");

        if (geometry.contains("dVoxel")) {
            const std::vector<double> dVoxel = geometry.at("dVoxel").get<std::vector<double>>();
            config.geometry.dVoxel.x = dVoxel[0];
            config.geometry.dVoxel.y = dVoxel[1];
        }
        else
            throw std::runtime_error("In json_config load_config(const std::string&) could not find dVoxelvariable in geometry.");

        if (geometry.contains("nDetector"))
            config.geometry.nDetector = geometry.at("nDetector").get<int>();
        else
            throw std::runtime_error("In json_config load_config(const std::string&) could not find nDetector in geometry");
        
        if (geometry.contains("dDetector"))
            config.geometry.dDetector = geometry.at("dDetector").get<double>();
        else
            throw std::runtime_error("In json_config load_config(const std::string&) could not find dDetector in geometry");

        if (geometry.contains("sDetector"))
            config.geometry.sDetector = geometry.at("sDetector").get<double>();
        else
            throw std::runtime_error("In json config load_config(const std::string&) could not find sDetector in geometry");

        if (geometry.contains("minAngle") && 
            geometry.contains("maxAngle") && 
            geometry.contains("nAngles")) {
            double minAngle = geometry.at("minAngle").get<double>();
            double maxAngle = geometry.at("maxAngle").get<double>();
            int nAngles = geometry.at("nAngles").get<int>();

            if (nAngles < 1)
                throw std::runtime_error("In json_config load_config(const std::string&) invalid angle configuration.");
            else if  (minAngle < 0 || maxAngle < 0 || minAngle > 360 || maxAngle > 360)
                throw std::runtime_error("In json_config load_config(const std::string&) rotation angles should be contained between 0 and 360");
            else if (maxAngle > minAngle)
                throw std::runtime_error("In json_config load_config(const std::string&) maxAngle > minAngle");
            else if (nAngles == 1 && minAngle != maxAngle)
                throw std::runtime_error("In json_config load_config(const std::string&) invalid angle configuration.");

        
            config.geometry.nAngles = nAngles;
            config.geometry.angles.resize(nAngles);
            if (nAngles == 1)
                config.geometry.angles[0] = minAngle;
            else {
                for (int i = 0; i < nAngles; i++) {
                    config.geometry.angles[i] = minAngle + (maxAngle - minAngle) / (nAngles-1) * i;
                }
            }
        }
        else
            throw std::runtime_error("In json config load_config(const std::string&) could not find rotation angle configuration in geometry");
    }

    // Done
    return config;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Main function
///////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {
    // Check input arguments
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " config.json\n";
        return 1;
    }

    // Load configuration
    json_config config = load_config(argv[1]);

    // Load image
    std::ifstream file(config.path_image.c_str(), std::ios::binary);
    if (!file)
        throw std::runtime_error("Could not open image file.");

    // Read the image
    opticalflow::Image image(config.geometry.nVoxel);
    file.read(reinterpret_cast<char*>(image.get_field()), image.get_size()*sizeof(double));

    if (!file)
        throw std::runtime_error("Could not read image file");

    // Allocate memory for projections
    const dim dimout(
        static_cast<std::size_t>(config.geometry.nDetector),
        static_cast<std::size_t>(config.geometry.nAngles)
    );
    opticalflow::Image projections(dimout);

    // Calculate the forward projection
    const auto start = std::chrono::steady_clock::now();
    fan_beam::Ax(projections, image, config.geometry);
    const auto end = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsed = end-start;
    std::cout << "Runtime (s): " << elapsed.count() << "\n";

    // Write to output file
    std::ofstream output(config.path_projections.c_str(), std::ios::binary);
    if (!output)
        throw std::runtime_error("Could not open output file.");
    output.write(reinterpret_cast<const char*>(projections.get_field()), projections.get_size()*sizeof(double));
    if (!output)
        throw std::runtime_error("Could not open output file.");
    
    // Done
    return 0;
}