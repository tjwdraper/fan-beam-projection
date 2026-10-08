#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <map>

#include "coord2d.hpp"
#include "json.hpp"
#include "Field.hpp"
#include "fan_beam.hpp"

///////////////////////////////////////////////////////////////////////////////////////////////////
// Read json configuration file and store in structure
///////////////////////////////////////////////////////////////////////////////////////////////////
struct json_config {
    // Path to files
    std::string path_image;
    std::string path_projections;

    double DSD;
    double DSO;
    
    dim nVoxel;
    vector2d sVoxel;
    vector2d dVoxel;

    int nDetector;
    double sDetector;
    double dDetector;
};

json_config load_config(const std::string& filename) {
    // Open .json configuration file
    std::ifstream file(filename);
    if (!file.is_open())
        throw std::runtime_error("Could not open json configuration fiel: " + filename);
    nlohmann::json json;
    file >> json;

    // Create json configuration structure
    json_config config;

    // Set the image and projection paths
    config.path_image = json.at("path_image").get<std::string>();
    config.path_projections = json.at("path_projections").get<std::string>();

    // Set fan beam geometry
    const auto& geometry = json.at("geometry");
    config.DSO = geometry.at("DSO").get<double>();
    config.DSD = geometry.at("DSD").get<double>();

    const auto nVoxel = geometry.at("nVoxel").get<std::array<int,2>>();
    config.nVoxel.x = nVoxel[0];
    config.nVoxel.y = nVoxel[1];

    const auto sVoxel = geometry.at("sVoxel").get<std::array<double,2>>();
    config.sVoxel.x = sVoxel[0];
    config.sVoxel.y = sVoxel[1];

    const auto dVoxel = geometry.at("dVoxel").get<std::array<double,2>>();
    config.dVoxel.x = dVoxel[0];
    config.dVoxel.y = dVoxel[1];

    // Detector
    config.nDetector = geometry.at("nDetector").get<int>();
    config.sDetector = geometry.at("sDetector").get<double>();
    config.dDetector = geometry.at("dDetector").get<double>();

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

}