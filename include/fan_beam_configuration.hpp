#ifndef _FAN_BEAM_CONFIGURATION_HPP_
#define _FAN_BEAM_CONFIGURATION_HPP_

#include <vector>
#include "coord2d.hpp"

class fan_beam_configuration {
    public:
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

        // Constructor
        fan_beam_configuration() = default;

        // Obtain location of fan-beam source and detector pixel under a rotating gantry.
        vector2d get_source_position(double angle) const {
            return rotate(vector2d(DSO, 0), angle);
        }

        vector2d get_detector_pixel_position(int u, double angle) const {
            if (u < 0 || u > nDetector-1)
                throw std::runtime_error("detector pixel index out of scope.");
            return rotate(vector2d(-(DSD-DSO), (u+0.5)*dDetector - sDetector/2), angle);
        }

        // Convert voxel index to physical coordinate (mm)
        vector2d voxel_to_mm(const vector2d& pos) const {
            double cx = static_cast<double>(nVoxel.x) / 2.0 - 0.5;
            double cy = static_cast<double>(nVoxel.y) / 2.0 - 0.5;
            return vector2d(
                (pos.x - cx) * dVoxel.x,
                (pos.y - cy) * dVoxel.y
            );
        }

        // Inverse: convert physical coordinate to image index.
        vector2d mm_to_voxel(const vector2d& pos) const {
            double cx = static_cast<double>(nVoxel.x) / 2.0 - 0.5;
            double cy = static_cast<double>(nVoxel.y) / 2.0 - 0.5;
            return vector2d(
                pos.x / dVoxel.x + cx,
                pos.y / dVoxel.y + cy
            );
        }

        // Compute distances between source and detector pixels
        std::vector<double> compute_ray_lengths() const {
            std::vector<double> raylengths(nDetector);

            const vector2d position_source = get_source_position(0);
            for (int u = 0; u < nDetector; u++) {
                const vector2d position_detector_pixel = get_detector_pixel_position(u, 0);
                const vector2d dir = position_detector_pixel - position_source;

                raylengths[u] = std::hypot(dir.x, dir.y);
            }

            return raylengths;
        }

        // Get physical location of the edges of the image
        double xmin() const {
            return -(static_cast<double>(nVoxel.x) / 2) * dVoxel.x;
        }

        double xmax() const {
            return (static_cast<double>(nVoxel.x) / 2) * dVoxel.x;
        }

        double ymin() const {
            return -(static_cast<double>(nVoxel.y) / 2) * dVoxel.y;
        }

        double ymax() const {
            return static_cast<double>(nVoxel.y) / 2 * dVoxel.y;
        }
};

#endif