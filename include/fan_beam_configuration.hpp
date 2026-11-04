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

    // Methods
    vector2d get_source_position(double angle) const {
        return rotate(vector2d(DSO, 0), angle);
    }

    vector2d get_detector_pixel_position(int u, double angle) const {
        if (u < 0 || u > nDetector-1)
            throw std::runtime_error("detector pixel index out of scope.");
        return rotate(vector2d(-(DSD-DSO), (u+0.5)*dDetector - sDetector/2), angle);
    }
};

#endif