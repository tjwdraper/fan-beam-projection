#ifndef _FAN_BEAM_HPP_
#define _FAN_BEAM_HPP_

#include "coord2d.hpp"
#include "Field.hpp"
#include "fan_beam_configuration.hpp"


#include <cmath>
#include <algorithm>
#include <limits>

enum class Dimension {X, Y};

inline int sgn(const double x) {
    return (x > 0) - (x < 0);
}

inline void clip_axis(double& tEnter, double& tExit, const double s, const double d, const double lo, const double hi) {
    double a = (lo - s) / d;
    double b = (hi - s) / d;

    if (a > b)
        std::swap(a,b);

    tEnter = std::max(tEnter, a);
    tExit = std::min(tExit, b);

    if (tEnter > tExit)
        throw std::runtime_error("In clip_axis(double&, double&, const double, const double, const double, const double): tEnter > tExit.");
}

namespace fan_beam {
    // Rotation of 2d position
    inline vector2d rotate(vector2d vec, double angle) {
        return vector2d(
            vec.x * cos(angle) - vec.y * sin(angle),
            vec.x * sin(angle) + vec.y * cos(angle)
        );
    }

    // Helper functions
    inline vector2d get_source_position(double angle, const fan_beam_configuration& config) {
        return rotate(vector2d(config.DSO,0), angle);
    }

    inline vector2d get_detector_pixel_position(int u, double angle, const fan_beam_configuration& config) {
        return rotate(vector2d(-(config.DSD-config.DSO),(u+0.5)*config.dDetector -config.sDetector/2.0), angle);
    }

    // Traverse a ray and call visit(x, y, segment_length) for every
    // voxel crossed. Assumes sVoxel = nVoxel * dVoxel.
    inline double accumulate_ray(vector2d source, vector2d dir, const opticalflow::Image& image, const fan_beam_configuration& config) {
        // Get pixel dimensions
        const int nx = static_cast<int>(config.nVoxel.x);
        const int ny = static_cast<int>(config.nVoxel.y);

        const double dVoxelX = config.dVoxel.x;
        const double dVoxelY = config.dVoxel.y;

        const double hx = 0.5 * config.sVoxel.x;
        const double hy = 0.5 * config.sVoxel.y;

        if (nx <= 0 || ny <= 0 || dVoxelX <= 0 || dVoxelY <= 0)
            return 0.0;

        const double rayLength = std::hypot(dir.x, dir.y);
        if (rayLength == 0.0)
            return 0.0;

        




        // Calculate the ray orientation. Does it travel along the positive or negative x/y-dimensions, and which dimension is dominant?
        Dimension dominantAxis = (std::abs(dir.x) / dVoxelX > std::abs(dir.y) / dVoxelY ? Dimension::X : Dimension::Y);
        int rayOrientationX = sgn(dir.x);
        int rayOrientationY = sgn(dir.y);

        // Compute the physical location when the ray enters the image volume
        const double xmin = -hx;
        const double xmax = hx;
        const double ymin = -hy;
        const double ymax = hy;

        double px_mm, py_mm;
        if (rayOrientationX == 1)
            px_mm = std::max(xmin, source.x);
        else if (rayOrientationX == -1)
            px_mm = std::min(xmax, source.x);
        else {
            if (source.x < xmin || source.x > xmax)
                throw std::runtime_error("If ray is perpendicular, it has to be contained in image dimensions to intersect.");
            px_mm = source.x;
        }

        if (rayOrientationY == 1)
            py_mm = std::max(ymin, source.y);
        else if (rayOrientationY == -1)
            py_mm = std::min(ymax, source.y);
        else {
            if (source.y < ymin || source.y > ymax)
                throw std::runtime_error("If ray is perpendicular, it has to be contained in image dimensions to intersect.");
            py_mm = source.y;
        }


        double px = (px_mm + hx) / dVoxelX;
        double py = (py_mm + hy) / dVoxelY;

        int dx = static_cast<int>(std::floor(px));
        int dy = static_cast<int>(std::floor(py));

        if (dx < 0 || dx > nx-1 || dy < 0 || dy > ny-1)
            throw std::runtime_error("In accumulate_ray, invalid entry location of ray in image volume");

        // Compute projection
        double p(0.0);

        // Iterate through x-dimension, interpolate over y-axis
        if (dominantAxis == Dimension::X) {
            while (dx > 0 && dx < nx-1 && dy > 0 && dy < ny-1) {
                const double f = py - dy;

                p += (1-f) * image.get_val(dx, dy) + f * image.get_val(dx, std::min(dy+1, ny-1));

                // Increment
                px += rayOrientationX;
                py += dir.y / dir.x;

                dx = static_cast<int>(std::floor(px));
                dy = static_cast<int>(std::floor(py));
            }
        }
        else if (dominantAxis == Dimension::Y) {
            while (dx > 0 && dx < nx-1 && dy > 0 && dy < ny-1) {
                const double f = px - dx;

                p += (1-f) * image.get_val(dx, dy) + f * image.get_val(std::min(dx+1, nx-1), dy);

                // Increment
                px += dir.x / dir.y;
                py += rayOrientationY;

                dx = static_cast<int>(std::floor(px));
                dy = static_cast<int>(std::floor(py));
            }
        }

        // Done
        return p;
    }

    template <typename Visit>
    inline void trace_ray(
        vector2d source, 
        vector2d dir, 
        double rayLength,
        const fan_beam_configuration& config,
        Visit visit
    ) {
        // Get physical location of image edge (in mm)
        const double cx = 0.5*config.sVoxel.x; //config.dVoxel.x * (0.5 * config.nVoxel.x);
        const double cy = 0.5*config.sVoxel.y; //config.dVoxel.y * (0.5 * config.nVoxel.y);
        const double xmin = -cx; const double xmax = cx;
        const double ymin = -cy; const double ymax = cy;

        // Calculate the values of t when the ray enters the image domain, and when it exits
        double tEnter = 0.0;
        double tExit = 1.0;

        clip_axis(tEnter, tExit, source.x, dir.x, xmin, xmax);
        clip_axis(tEnter, tExit, source.y, dir.y, ymin, ymax);

        if (tEnter > tExit)
            return;

        // Get the orientation of the ray: is it moving in the positive/negative x/y dimensions
        double beamOrientationX = sgn(dir.x);
        double beamOrientationY = sgn(dir.y);

        // Calculate the position where the ray enters the image. Use std::nextafter to infinitesimally move into the image domain.
        double px = std::nextafter(source.x + tEnter * dir.x, source.x + (tEnter + 1.0) * dir.x);
        double py = std::nextafter(source.y + tEnter * dir.y, source.y + (tEnter + 1.0) * dir.y);
        
        // Convert to voxel coordinates, clip to the image boundary
        int dx = static_cast<int>(std::floor((px + cx) / config.dVoxel.x));
        int dy = static_cast<int>(std::floor((py + cy) / config.dVoxel.y));

        dx = std::min(static_cast<int>(config.nVoxel.x-1), std::max(0, dx));
        dy = std::min(static_cast<int>(config.nVoxel.y-1), std::max(0, dy));

        // Compute time it takes to traverse of voxel in the x- and y-direction. If the ray is parallel to the cooridinate axes 
        // (which is the case if beamOrientationXY == 0), it takes infinite time to traverse between voxel boundaries.
        const double xdt = beamOrientationX == 0 ? std::numeric_limits<double>::infinity() : config.dVoxel.x / std::abs(dir.x);
        const double ydt = beamOrientationY == 0 ? std::numeric_limits<double>::infinity() : config.dVoxel.y / std::abs(dir.y);

        // Compute the time the ray passes through the next boundary as seen from the first voxel entry
        double tNextX = std::numeric_limits<double>::infinity();
        double tNextY = std::numeric_limits<double>::infinity();

        if (beamOrientationX != 0) {
            const double boundary = xmin + (beamOrientationX > 0 ? dx + 1 : dx) * config.dVoxel.x;
            tNextX = (boundary - source.x) / dir.x;
        }
        if (beamOrientationY != 0) {
            const double boundary = ymin + (beamOrientationY > 0 ? dy + 1 : dy) * config.dVoxel.y;
            tNextY = (boundary - source.y) / (dir.y);
        }

        // Move through the participating voxels in the image domain
        double tCur = tEnter;
        while(tCur < tExit &&
            dx >= 0 && dx <= config.nVoxel.x-1 &&
            dy >= 0 && dy <= config.nVoxel.y-1) {

            // Compute the time of the next event
            const double tNext = std::min({tNextX, tNextY, tExit});

            // Compute contribution of current voxel
            if (tNext > tCur)
                visit(dx, dy, (tNext - tCur) * rayLength);

            // Update variables
            if (tNextX == tNext) {
                dx += beamOrientationX;
                tNextX += xdt;
            }

            if (tNextY == tNext) {
                dy += beamOrientationY;
                tNextY += ydt;
            }

            // Update tCur
            tCur = tNext;
        }
    }

    // Forward projection
    inline void Ax(opticalflow::Image& projections,
        const opticalflow::Image& image,
        const fan_beam_configuration& config
    ) {
        // Iterate over rotation angles
        for (std::size_t a = 0; a < config.angles.size(); ++a) {
            // Get the current rotation angle
            const double angle = config.angles[a];

            // Calculate position of the source
            vector2d position_source = fan_beam::get_source_position(angle, config);

            for (std::size_t u = 0; u < config.nDetector; ++u) {
                // Get the position of the voxel in mm and the direction vector between the source and current pixel:
                vector2d position_detector_pixel = fan_beam::get_detector_pixel_position(u, angle, config);
                vector2d dir = position_detector_pixel - position_source;

                // Get the length of the ray
                const double DSP = std::hypot(dir.x, dir.y);

                // Calculate the contribution of the ray
                // double sum = 0.0;                
                // fan_beam::trace_ray(position_source, dir, DSP, config, 
                //     [&](int x, int y, double length) {
                //         sum += length * image.get_val(x, y);
                //     }
                // );

                double sum = fan_beam::accumulate_ray(position_source, dir, image, config);

                projections.set_val(sum, u, a);
            }
        }
    }

    inline void Atb(opticalflow::Image& backproj,
        const opticalflow::Image& projections,
        const fan_beam_configuration& config
    ) {
        // Iterate over rotation angles
        for (std::size_t a = 0; a < config.angles.size(); ++a) {
            // Get the current rotation angle
            const double angle = config.angles[a];

            // Calculate position of the 
            vector2d position_source = fan_beam::get_source_position(angle, config);

            for (std::size_t u = 0; u < config.nDetector; ++u) {
                // Get the position of the voxel in mm and the direction vector between the source and current pixel:
                vector2d position_detector_pixel = fan_beam::get_detector_pixel_position(u, angle, config);
                vector2d dir = position_detector_pixel - position_source;

                // Get the length of the ray
                const double DSP = std::hypot(dir.x, dir.y);
                const double p = projections.get_val(u,a);

                // Calculate the contribution of the ray
                fan_beam::trace_ray(position_source, dir, DSP, config, 
                    [&](int x, int y, double length) {
                        backproj.set_val(backproj.get_val(x,y) + length*p, x, y);
                    }
                );                
            }

        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Primal-dual methods
    ///////////////////////////////////////////////////////////////////////////////////////////////
    
    void update_p(opticalflow::Image& px, opticalflow::Image& py,
                  const opticalflow::Image& ubar, 
                  const opticalflow::Image& vbarx, const opticalflow::Image& vbary,
                  double sigma) {
        const dim dimin = px.get_dimensions();

        for (std::size_t j = 0; j < dimin.y; ++j) {
            for (std::size_t i = 0; i < dimin.x; ++i) {
                px.set_val(px.get_val(i,j) + sigma * (opticalflow::gradients::dx(ubar, i, j) - vbarx.get_val(i,j)),i,j);
                py.set_val(py.get_val(i,j) + sigma * (opticalflow::gradients::dy(ubar, i, j) - vbary.get_val(i,j)),i,j);
            }
        }
    }

    void update_q(opticalflow::Image& qxx, opticalflow::Image& qyy, opticalflow::Image& qxy,
                  const opticalflow::Image& vbarx, const opticalflow::Image& vbary,
                  double sigma) {
        const dim dimin = qxx.get_dimensions();

        for (std::size_t j = 0; j < dimin.y; ++j) {
            for (std::size_t i = 0; i < dimin.x; ++i) {
                const double Evbarxx = opticalflow::gradients::dx(vbarx, i, j);
                const double Evbaryy = opticalflow::gradients::dy(vbary, i, j);
                const double Evbarxy = 0.5 * (opticalflow::gradients::dx(vbary,i,j) + 
                                              opticalflow::gradients::dy(vbarx,i,j));

                qxx.set_val(qxx.get_val(i,j) + sigma * Evbarxx,i,j);
                qyy.set_val(qyy.get_val(i,j) + sigma * Evbaryy,i,j);
                qxy.set_val(qxy.get_val(i,j) + sigma * Evbarxy,i,j);
            }
        }
    }

    void update_r(opticalflow::Image& r, const opticalflow::Image& ubar, const opticalflow::Image& p, double sigma,
        const fan_beam_configuration& config) {
        const int nAngles = config.angles.size();

        opticalflow::Image Kr(dim(config.nDetector, config.nAngles)); Kr.fill(0.0);
        fan_beam::Ax(Kr, ubar, config);
        Kr/= 300;

        for (std::size_t a = 0; a < config.nAngles; ++a) {
            for (std::size_t u = 0; u < config.nDetector; ++u) {
                r.set_val(r.get_val(u,a) + sigma * (Kr.get_val(u,a) - p.get_val(u,a)), u,a);
            }
        }
    }

    void update_u(opticalflow::Image& u,
                  const opticalflow::Image& px, const opticalflow::Image& py, const opticalflow::Image r,
                  double tau,
                  const fan_beam_configuration& config) {
        const dim dimin = u.get_dimensions();        

        opticalflow::Image Ktr(dimin); Ktr.fill(0.0);
        fan_beam::Atb(Ktr, r, config);
        Ktr /= 300;

        for (std::size_t j = 0; j < dimin.y; ++j) {
            for (std::size_t i = 0; i < dimin.x; ++i) {
                const double divp = opticalflow::gradients::adjoint_dx(px,i,j) + opticalflow::gradients::adjoint_dy(py,i,j);
                u.set_val(u.get_val(i,j) + tau*(divp - Ktr.get_val(i,j)),i,j);
            }
        }

    }

    void update_v(opticalflow::Image& vx, opticalflow::Image& vy,
                  const opticalflow::Image& px, const opticalflow::Image& py,
                  const opticalflow::Image& qxx, const opticalflow::Image& qyy, const opticalflow::Image& qxy,
                  double tau) {
        const dim dimin = vx.get_dimensions();

        for (std::size_t j = 0; j < dimin.y; ++j) {
            for (std::size_t i = 0; i < dimin.x; ++i) {
                const double divqx = opticalflow::gradients::adjoint_dx(qxx,i,j) + opticalflow::gradients::adjoint_dy(qxy,i,j);
                const double divqy = opticalflow::gradients::adjoint_dy(qyy,i,j) + opticalflow::gradients::adjoint_dx(qxy,i,j);

                vx.set_val(vx.get_val(i,j) + tau * (px.get_val(i,j) + divqx), i,j);
                vy.set_val(vy.get_val(i,j) + tau * (py.get_val(i,j) + divqy), i,j);
            }
        }

    }

    void proj_p(opticalflow::Image& px, opticalflow::Image& py, double alpha) {
        const dim dimin = px.get_dimensions();

        for (std::size_t j = 0; j < dimin.y; ++j) {
            for (std::size_t i = 0; i < dimin.x; ++i) {
                const double pxv = px.get_val(i,j);
                const double pyv = py.get_val(i,j);

                const double normp = std::sqrt(pxv*pxv + pyv*pyv) / alpha;

                if (normp > 1.0) {
                    px.set_val(px.get_val(i,j)/normp, i,j);
                    py.set_val(py.get_val(i,j)/normp, i,j);
                }
            }
        }
    }

    void proj_q(opticalflow::Image& qxx, opticalflow::Image& qyy, opticalflow::Image& qxy, double alpha) {
        const dim dimin = qxx.get_dimensions();

        for (std::size_t j = 0; j < dimin.y; ++j) {
            for (std::size_t i = 0; i < dimin.x; ++i) {
                const double qxxv = qxx.get_val(i,j);
                const double qyyv = qyy.get_val(i,j);
                const double qxyv = qxy.get_val(i,j);

                const double normq = std::sqrt(qxxv*qxxv + qyyv*qyyv + 2*qxyv*qxyv) / alpha;

                if (normq > 1.0) {
                    qxx.set_val(qxx.get_val(i,j)/normq, i,j);
                    qyy.set_val(qyy.get_val(i,j)/normq, i,j);
                    qxy.set_val(qxy.get_val(i,j)/normq, i,j);
                }
            }
        }
    }

    void prox_r(opticalflow::Image& r, double sigma, double lambda) {
        const dim dimin = r.get_dimensions();

        for (std::size_t a = 0; a < dimin.y; ++a) {
            for (std::size_t u = 0; u < dimin.x; ++u) {
                r.set_val(r.get_val(u,a) / (1.0 + lambda * sigma), u, a);
            }
        }
    }

    // The primal-dual TGV-denoising algorithm
    opticalflow::Image tgv_reconstruct(const opticalflow::Image& p, 
        double tau, double sigma, double lambda, double alpha0, double alpha1, int niter, double convergence, 
        const fan_beam_configuration& config) {
        
        const int nAngles = config.angles.size();

        // Initialize Primal-Dual variables
        opticalflow::Image u(config.nVoxel); u.fill(0.0);
        opticalflow::Image ubar(config.nVoxel); ubar.fill(0.0);

        opticalflow::Image vx(config.nVoxel); vx.fill(0.0);
        opticalflow::Image vy(config.nVoxel); vy.fill(0.0);

        opticalflow::Image vbarx(config.nVoxel); vbarx.fill(0.0);
        opticalflow::Image vbary(config.nVoxel); vbary.fill(0.0);

        opticalflow::Image px(config.nVoxel); px.fill(0.0);
        opticalflow::Image py(config.nVoxel); py.fill(0.0);

        opticalflow::Image qxx(config.nVoxel); qxx.fill(0.0);
        opticalflow::Image qyy(config.nVoxel); qyy.fill(0.0);
        opticalflow::Image qxy(config.nVoxel); qxy.fill(0.0);

        opticalflow::Image r(dim(config.nDetector, config.nAngles)); r.fill(0.0);

        // Tracking variables
        opticalflow::Image uold(config.nVoxel); uold.fill(0.0);
        opticalflow::Image vxold(config.nVoxel); vxold.fill(0.0);
        opticalflow::Image vyold(config.nVoxel); vyold.fill(0.0);    

        // Primal-dual iterations
        for (int iter = 0; iter < niter; ++iter) {
            // Update dual variables
            fan_beam::update_p(px,py,ubar,vbarx,vbary,sigma);
            fan_beam::update_q(qxx, qyy, qxy, vbarx, vbary, sigma);
            fan_beam::update_r(r, ubar, p, sigma, config);

            fan_beam::proj_p(px,py,alpha1);
            fan_beam::proj_q(qxx,qyy,qxy,alpha0);
            fan_beam::prox_r(r, sigma, lambda);

            // Track u
            uold = u;

            // Proximal operator
            fan_beam::update_u(u, px, py, r, tau, config);

            // Update ubar
            ubar = 2*u - uold;

            // Track v
            vxold = vx;
            vyold = vy;

            // Update v
            fan_beam::update_v(vx, vy, px, py, qxx, qyy, qxy, tau);

            // Update vbar
            vbarx = 2*vx - vxold;
            vbary = 2*vy - vyold;

            // Check for convergence
            double relchange = opticalflow::image::norm(u-uold)/opticalflow::image::norm(u);
            if (relchange < convergence)
                break;

            // Update some norms:
            if (iter % 50 == 0) {
                std::cout << "iter " << iter 
                    << "\t|u_new-u|/|u_new| = " << relchange 
        //             // << "\t1/(2*lambda)*|u-f|^2 = " << opticalflow::image::normsq(u-f) / (2*lambda)
        //             // << "\talpha_1 * |Du - v|_1 = " << denoise::tv_residual_norm(u,vx,vy) * alpha1
        //             // << "\talpha_0 * |Ev|_1 = " << denoise::l1_E_norm(vx,vy) * alpha0
                    << std::endl;
            }
        }

        return u;
    }
}

#endif