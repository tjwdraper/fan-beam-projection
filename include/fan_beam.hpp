#ifndef _FAN_BEAM_HPP_
#define _FAN_BEAM_HPP_

#include "coord2d.hpp"
#include "Field.hpp"



#include <cmath>
#include <algorithm>
#include <limits>

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

    // Traverse a ray and call visit(x, y, segment_length) for every
    // voxel crossed. Assumes sVoxel = nVoxel * dVoxel.
    template <typename Visit>
    inline void trace_ray(
        vector2d source,
        vector2d dir,
        dim nVoxel,
        vector2d sVoxel,
        vector2d dVoxel,
        double rayLength,
        Visit visit
    ) {
        const double xmin = -0.5 * sVoxel.x;
        const double ymin = -0.5 * sVoxel.y;
        const double xmax = xmin + nVoxel.x * dVoxel.x;
        const double ymax = ymin + nVoxel.y * dVoxel.y;

        // Clip ray parameter t to [0, 1] and the volume bounding box.
        double tEnter = 0.0;
        double tExit  = 1.0;

        auto clip_axis = [&](double s, double d,
                            double lo, double hi) -> bool {
            if (std::abs(d) < 1e-14) {
                return s >= lo && s <= hi;
            }

            double a = (lo - s) / d;
            double b = (hi - s) / d;
            if (a > b) std::swap(a, b);

            tEnter = std::max(tEnter, a);
            tExit  = std::min(tExit, b);

            return tEnter < tExit;
        };

        if (!clip_axis(source.x, dir.x, xmin, xmax)) return;
        if (!clip_axis(source.y, dir.y, ymin, ymax)) return;

        if (!(tEnter < tExit)) return;

        // Move infinitesimally into the volume to select the correct
        // voxel when the entry point lies exactly on a voxel boundary.
        const double px = std::nextafter(
            source.x + tEnter * dir.x,
            source.x + (tEnter + 1.0) * dir.x
        );
        const double py = std::nextafter(
            source.y + tEnter * dir.y,
            source.y + (tEnter + 1.0) * dir.y
        );

        int ix = static_cast<int>(std::floor((px - xmin) / dVoxel.x));
        int iy = static_cast<int>(std::floor((py - ymin) / dVoxel.y));

        ix = std::max(0, std::min(ix, static_cast<int>(nVoxel.x) - 1));
        iy = std::max(0, std::min(iy, static_cast<int>(nVoxel.y) - 1));

        const int stepX = (dir.x > 0.0) - (dir.x < 0.0);
        const int stepY = (dir.y > 0.0) - (dir.y < 0.0);

        const double inf = std::numeric_limits<double>::infinity();

        const double tDeltaX =
            stepX == 0 ? inf : dVoxel.x / std::abs(dir.x);
        const double tDeltaY =
            stepY == 0 ? inf : dVoxel.y / std::abs(dir.y);

        double tMaxX = inf;
        double tMaxY = inf;

        if (stepX != 0) {
            const double boundary = xmin +
                (stepX > 0 ? ix + 1 : ix) * dVoxel.x;
            tMaxX = (boundary - source.x) / dir.x;
        }

        if (stepY != 0) {
            const double boundary = ymin +
                (stepY > 0 ? iy + 1 : iy) * dVoxel.y;
            tMaxY = (boundary - source.y) / dir.y;
        }

        double t = tEnter;

        while (t < tExit && ix >= 0 && ix < nVoxel.x &&
            iy >= 0 && iy < nVoxel.y) {

            const double tNext = std::min({tMaxX, tMaxY, tExit});

            if (tNext > t) {
                visit(ix, iy, (tNext - t) * rayLength);
            }

            // Advance across every boundary reached at this parameter.
            // Advancing both axes handles a ray passing through a corner.
            if (tMaxX <= tNext) {
                ix += stepX;
                tMaxX += tDeltaX;
            }

            if (tMaxY <= tNext) {
                iy += stepY;
                tMaxY += tDeltaY;
            }

            t = tNext;
        }
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

            // Calculate position of the source
            vector2d position_source = get_source_position(angle, DSO);

            for (std::size_t u = 0; u < nDetector; ++u) {
                // Get the position of the voxel in mm and the direction vector between the source and current pixel:
                vector2d position_detector_pixel = get_detector_pixel_position(u, angle, DSO, DSD, dDetector, sDetector);
                vector2d dir = position_detector_pixel - position_source;

                // Get the length of the ray
                const double DSP = std::hypot(dir.x, dir.y);

                // Calculate the contribution of the ray
                double sum = 0.0;
                trace_ray(position_source, dir, nVoxel, sVoxel, dVoxel, DSP, 
                    [&](int x, int y, double length) {
                        sum += length*image.get_val(x,y);
                    }
                );

                projections.set_val(sum, u, a);

                // // Calculate the ray intersecting with voxel boundaries
                // std::vector<double> t = get_events(position_source, dir, nVoxel, sVoxel, dVoxel);

                // // Iterate over time intervals
                // double p = accumulate_ray(t, image, position_source, dir, sVoxel, dVoxel, DSP);
                // projections.set_val(p,u,a);
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

                // Get the length of the ray
                const double DSP = std::hypot(dir.x, dir.y);
                const double p = projections.get_val(u,a);

                // Calculate the contribution of the ray
                trace_ray(position_source, dir, nVoxel, sVoxel, dVoxel, DSP, 
                    [&](int x, int y, double length) {
                        backproj.set_val(backproj.get_val(x,y) + length*p, x, y);
                    }
                );

                // // Calculate the ray intersecting with voxel boundaries
                // std::vector<double> events = get_events(position_source, dir, nVoxel, sVoxel, dVoxel);

                // // Calculate distance between source and pixel
                // double DSP = std::sqrt(dir.x*dir.x + dir.y*dir.y);

                // // Add contribution of current projection to back-projection
                // fan_beam::adjoint_accumulate_ray(backproj, events, projections.get_val(u,a), position_source, dir, sVoxel, dVoxel, DSP);
                
            }

        }
    }
}

#endif