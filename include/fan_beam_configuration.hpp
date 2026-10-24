#ifndef _FAN_BEAM_CONFIGURATION_HPP_
#define _FAN_BEAM_CONFIGURATION_HPP_

#include <vector>
#include "coord2d.hpp"

struct fan_beam_configuration {
    double DSO; // Distance source to origin
    double DSD; // Distance source to detector

    dim nVoxel; // Number of voxels
    vector2d sVoxel; // Size of image in mm
    vector2d dVoxel; // Voxel size, i.e. dVoxel = sVoxel / nVoxel
    
    int nDetector; // Number of pixels on detector panel
    double sDetector; // Size of detector in mm
    double dDetector; // Pixel size of detector, i.e. dDetector = sDetector / nDetector

    std::vector<double> angles; // gantry angles or fan-beam acquisition angles.
    int nAngles;
};

#endif