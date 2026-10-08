#ifndef _FAN_BEAM_HPP_
#define _FAN_BEAM_HPP_

#include "coord2d.hpp"
#include "Field.hpp"

namespace fan_beam {
    // Rotation of 2d position
    inline vector2d rotate(vector2d vec, double angle) {
        return vector2d(
            vec.x * cos(angle) - vec.y * sin(angle),
            vec.x * sin(angle) + vec.y * cos(angle)
        );
    }

    // Helper functions
    inline vector2d get_source_position(double angle, double DSO) {
        return rotate(vector2d(DSO,0), angle);
    }

    inline vector2d get_detector_pixel_position(int u, double angle, double DSO, double DSD, double dDetector, double sDetector) {
        return rotate(vector2d(-(DSD-DSO),(u+0.5)*dDetector -sDetector/2.0), angle);
    }

    inline std::vector<double> get_events(vector2d position_source, vector2d dir, dim nVoxel, vector2d sVoxel, vector2d dVoxel) {
        // Allocate memory to check 
        std::vector<double> tx(nVoxel.x+1); // Nx voxels -> Nx+1 boundaries
        std::vector<double> ty(nVoxel.y+1); // Ny voxels -> Ny+1 boundaries

        for (std::size_t i = 0; i < nVoxel.x+1; ++i)
            tx[i] = (i*dVoxel.x - sVoxel.x/2 - position_source.x) / dir.x;

        for (std::size_t j = 0; j < nVoxel.y+1; ++j)
            ty[j] = (j*dVoxel.y - sVoxel.y/2 - position_source.y) / dir.y;

        // Combine unique values of tx, ty
        std::vector<double> t;
        for (double v : tx) {
            if (v >= 0 && v <= 1)
                t.push_back(v);
        }

        for (double v : ty) {
            if (v >= 0 && v <= 1)
                t.push_back(v);
        }

        std::sort(t.begin(), t.end());
        t.erase(std::unique(t.begin(), t.end()), t.end());

        return t;
    }

    inline double accumulate_ray(std::vector<double>& events, const opticalflow::Image& image, vector2d position_source, vector2d dir, vector2d sVoxel, vector2d dVoxel, double DSP) {
        double p = 0.0;
        const dim dimin = image.get_dimensions();
        for (std::size_t i = 0; i < events.size()-1; ++i) {
            double dt = events[i+1] - events[i];

            double mid_point = 0.5 * (events[i] + events[i+1]);
            int x = static_cast<int>(std::floor((position_source.x + mid_point * dir.x + sVoxel.x/2.0)/dVoxel.x));
            int y = static_cast<int>(std::floor((position_source.y + mid_point * dir.y + sVoxel.y/2.0)/dVoxel.y));

            if (x < 0 || x >= dimin.x || y < 0 || y >= dimin.y)
                continue;

            p += dt * DSP * image.get_val(x,y);
        }
        return p;
    }

    inline void adjoint_accumulate_ray(opticalflow::Image& backproj, const std::vector<double>& events, const double p, vector2d position_source, vector2d dir, vector2d sVoxel, vector2d dVoxel, double DSP) {
        const dim dimin = backproj.get_dimensions();
        for (std::size_t i = 0; i < events.size()-1; ++i) {
            double dt = events[i+1] - events[i];

            double mid_point = 0.5 * (events[i] + events[i+1]);
            int x = static_cast<int>(std::floor((position_source.x + mid_point * dir.x + sVoxel.x/2.0)/dVoxel.x));
            int y = static_cast<int>(std::floor((position_source.y + mid_point * dir.y + sVoxel.y/2.0)/dVoxel.y));

            if (x < 0 || x >= dimin.x || y < 0 || y >= dimin.y)
                continue;

            backproj.set_val(backproj.get_val(x,y) + dt * DSP * p,x,y);
        }
    }

    // Forward projection
    inline void Ax(opticalflow::Image& projections,
        const opticalflow::Image& image,
        double DSO, double DSD,
        dim nVoxel, vector2d sVoxel, vector2d dVoxel,
        int nDetector, double sDetector, double dDetector,
        const std::vector<double>& angles
    ) {
        // Iterate over rotation angles
        for (std::size_t a = 0; a < angles.size(); ++a) {
            // Get the current rotation angle
            const double angle = angles[a];

            // Calculate position of the 
            vector2d position_source = get_source_position(angle, DSO);

            for (std::size_t u = 0; u < nDetector; ++u) {
                // Get the position of the voxel in mm and the direction vector between the source and current pixel:
                vector2d position_detector_pixel = get_detector_pixel_position(u, angle, DSO, DSD, dDetector, sDetector);
                vector2d dir = position_detector_pixel - position_source;

                // Calculate the ray intersecting with voxel boundaries
                std::vector<double> t = get_events(position_source, dir, nVoxel, sVoxel, dVoxel);

                // Calculate distance between source and pixel
                double DSP = std::sqrt(dir.x*dir.x + dir.y*dir.y);

                // Iterate over time intervals
                double p = accumulate_ray(t, image, position_source, dir, sVoxel, dVoxel, DSP);
                projections.set_val(p,u,a);
            }
        }
    }

    inline void Atb(opticalflow::Image& backproj,
        const opticalflow::Image& projections,
        double DSO, double DSD,
        dim nVoxel, vector2d sVoxel, vector2d dVoxel,
        int nDetector, double sDetector, double dDetector,
        const std::vector<double>& angles
    ) {
        // Iterate over rotation angles
        for (std::size_t a = 0; a < angles.size(); ++a) {
            // Get the current rotation angle
            const double angle = angles[a];

            // Calculate position of the 
            vector2d position_source = get_source_position(angle, DSO);

            for (std::size_t u = 0; u < nDetector; ++u) {
                // Get the position of the voxel in mm and the direction vector between the source and current pixel:
                vector2d position_detector_pixel = get_detector_pixel_position(u, angle, DSO, DSD, dDetector, sDetector);
                vector2d dir = position_detector_pixel - position_source;

                // Calculate the ray intersecting with voxel boundaries
                std::vector<double> events = get_events(position_source, dir, nVoxel, sVoxel, dVoxel);

                // Calculate distance between source and pixel
                double DSP = std::sqrt(dir.x*dir.x + dir.y*dir.y);

                // Add contribution of current projection to back-projection
                fan_beam::adjoint_accumulate_ray(backproj, events, projections.get_val(u,a), position_source, dir, sVoxel, dVoxel, DSP);
                
            }

        }
    }
}

#endif