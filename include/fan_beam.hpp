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

inline bool clip_axis(double& tEnter, double& tExit, const double s, const double d, const double lo, const double hi) {
    if (std::abs(d) < 1e-12)
        return s >= lo && s <= hi;

    double a = (lo - s) / d;
    double b = (hi - s) / d;

    if (a > b)
        std::swap(a,b);

    tEnter = std::max(tEnter, a);
    tExit = std::min(tExit, b);

    return tEnter < tExit;
}

namespace fan_beam {
    // Rotation of 2d position
    

    

    // Traverse a ray and call visit(x, y, segment_length) for every
    // voxel crossed. Assumes sVoxel = nVoxel * dVoxel.
    inline double accumulate_ray(vector2d source, vector2d dir, const opticalflow::Image& image, const fan_beam_configuration& config) {
        // Get pixel dimensions
        const int nx = static_cast<int>(config.nVoxel.x);
        const int ny = static_cast<int>(config.nVoxel.y);

        // const double cx = static_cast<double>(nx)/2 - 0.5;
        // const double cy = static_cast<double>(ny)/2 - 0.5;

        const double dVoxelX = config.dVoxel.x;
        const double dVoxelY = config.dVoxel.y;

        // const double xmin = -(cx + 0.5) * dVoxelX;
        // const double xmax = (cx + 0.5) * dVoxelX;
        // const double ymin = -(cy + 0.5) * dVoxelY;
        // const double ymax = (cy + 0.5) * dVoxelY;

        if (nx <= 0 || ny <= 0 || dVoxelX <= 0 || dVoxelY <= 0)
            return 0.0;

        const double rayLength = std::hypot(dir.x, dir.y);
        if (rayLength == 0.0)
            return 0.0;

        // Calculate the ray orientation. Does it travel along the positive or negative x/y-dimensions, and which dimension is dominant?
        Dimension dominantAxis = (std::abs(dir.x) / dVoxelX > std::abs(dir.y) / dVoxelY ? Dimension::X : Dimension::Y);
        int rayOrientationX = sgn(dir.x);
        int rayOrientationY = sgn(dir.y);

        if (rayOrientationX == 0 && rayOrientationY == 0)
            throw std::runtime_error("Both ray orientations cannot be zero simultaneously.");

        // Compute intersection with the image rectangle.
        double tEnter = 0.0;
        double tExit  = 1.0;

        if (!clip_axis(tEnter, tExit, source.x, dir.x, config.xmin(), config.xmax()) ||
            !clip_axis(tEnter, tExit, source.y, dir.y, config.ymin(), config.ymax()) ||
            tEnter < 0 || tEnter > 1 ||
            tExit < 0 || tExit > 1 || 
            tEnter >= tExit)
        {
            return 0.0; // No positive-length intersection.
        }

        double p = 0.0;

        if (dominantAxis == Dimension::X)
        {
            if (std::abs(dir.x) < 1e-12)
                return 0.0;

            const double dt = dVoxelX / std::abs(dir.x);

            while (tEnter < tExit)
            {
                // Current integration interval.
                const double t0 = tEnter;
                const double t1 = std::min(t0 + dt, tExit);

                // Sample at the midpoint of the interval.
                const double tm = 0.5 * (t0 + t1);

                // const double x = source.x + tm * dir.x;
                // const double y = source.y + tm * dir.y;
                const vector2d pos_mm = source + tm * dir;

                // Continuous pixel-center coordinates.
                vector2d pos = config.mm_to_voxel(pos_mm);
                // const double px = x / dVoxelX + cx;
                // double py = y / dVoxelY + cy;

                // Major-axis pixel index.
                int ix = static_cast<int>(std::floor(pos.x + 0.5));
                ix = std::clamp(ix, 0, nx - 1);

                // Clamp minor coordinate to image edges.
                pos.y = std::clamp(pos.y, 0.0, static_cast<double>(ny - 1));

                const int iy0 = static_cast<int>(std::floor(pos.y));
                const int iy1 = std::min(iy0 + 1, ny - 1);

                const double fy = pos.y - iy0;

                const double value =
                    (1.0 - fy) * image.get_val(ix, iy0) +
                    fy         * image.get_val(ix, iy1);

                // Physical length represented by this sample.
                const double ds = (t1 - t0) * rayLength;

                p += value * ds;

                tEnter = t1;
            }
        }
        else
        {
            if (std::abs(dir.y) < 1e-12)
                return 0.0;

            const double dt = dVoxelY / std::abs(dir.y);

            while (tEnter < tExit)
            {
                const double t0 = tEnter;
                const double t1 = std::min(t0 + dt, tExit);

                const double tm = 0.5 * (t0 + t1);

                // const double x = source.x + tm * dir.x;
                // const double y = source.y + tm * dir.y;
                const vector2d pos_mm = source + tm * dir;

                // const double px = x / dVoxelX + cx;
                // const double py = y / dVoxelY + cy;
                vector2d pos = config.mm_to_voxel(pos_mm);

                int iy = static_cast<int>(std::floor(pos.y + 0.5));
                iy = std::clamp(iy, 0, ny - 1);

                double pxClamped = std::clamp(pos.x, 0.0, static_cast<double>(nx - 1));

                const int ix0 = static_cast<int>(std::floor(pxClamped));
                const int ix1 = std::min(ix0 + 1, nx - 1);

                const double fx = pxClamped - ix0;

                const double value =
                    (1.0 - fx) * image.get_val(ix0, iy) +
                    fx         * image.get_val(ix1, iy);

                const double ds = (t1 - t0) * rayLength;

                p += value * ds;

                tEnter = t1;
            }
        }

        return p;
    }

    inline void accumulate_ray_adjoint(opticalflow::Image& backproj, vector2d source, vector2d dir, const double p, const fan_beam_configuration& config) {
        // Get pixel dimensions
        const int nx = static_cast<int>(config.nVoxel.x);
        const int ny = static_cast<int>(config.nVoxel.y);

        const double cx = static_cast<double>(nx)/2 - 0.5;
        const double cy = static_cast<double>(ny)/2 - 0.5;

        const double dVoxelX = config.dVoxel.x;
        const double dVoxelY = config.dVoxel.y;

        const double xmin = -(cx + 0.5) * dVoxelX;
        const double xmax = (cx + 0.5) * dVoxelX;
        const double ymin = -(cy + 0.5) * dVoxelY;
        const double ymax = (cy + 0.5) * dVoxelY;

        if (nx <= 0 || ny <= 0 || dVoxelX <= 0 || dVoxelY <= 0)
            return;

        const double rayLength = std::hypot(dir.x, dir.y);
        if (rayLength == 0.0)
            return;
            
        // Calculate the ray orientation. Does it travel along the positive or negative x/y-dimensions, and which dimension is dominant?
        Dimension dominantAxis = (std::abs(dir.x) / dVoxelX > std::abs(dir.y) / dVoxelY ? Dimension::X : Dimension::Y);
        int rayOrientationX = sgn(dir.x);
        int rayOrientationY = sgn(dir.y);

        if (rayOrientationX == 0 && rayOrientationY == 0)
            throw std::runtime_error("Both ray orientations cannot be zero simultaneously.");

        // Compute intersection with the image rectangle.
        double tEnter = 0.0;
        double tExit  = 1.0;

        if (!clip_axis(tEnter, tExit, source.x, dir.x, xmin, xmax) ||
            !clip_axis(tEnter, tExit, source.y, dir.y, ymin, ymax) ||
            tEnter < 0 || tEnter > 1 ||
            tExit < 0 || tExit > 1 || 
            tEnter >= tExit)
        {
            return; // No positive-length intersection.
        }

        if (dominantAxis == Dimension::X)
        {
            if (std::abs(dir.x) < 1e-12)
                return;

            const double dt = dVoxelX / std::abs(dir.x);

            while (tEnter < tExit)
            {
                // Current integration interval.
                const double t0 = tEnter;
                const double t1 = std::min(t0 + dt, tExit);

                // Sample at the midpoint of the interval.
                const double tm = 0.5 * (t0 + t1);

                const double x = source.x + tm * dir.x;
                const double y = source.y + tm * dir.y;

                // Continuous pixel-center coordinates.
                const double px = x / dVoxelX + cx;
                double py = y / dVoxelY + cy;

                // Major-axis pixel index.
                int ix = static_cast<int>(std::floor(px + 0.5));
                ix = std::clamp(ix, 0, nx - 1);

                // Clamp minor coordinate to image edges.
                py = std::clamp(py, 0.0, static_cast<double>(ny - 1));

                const int iy0 = static_cast<int>(std::floor(py));
                const int iy1 = std::min(iy0 + 1, ny - 1);

                const double fy = py - iy0;

                // Physical length represented by this sample.
                const double ds = (t1 - t0) * rayLength;

                backproj.set_val(backproj.get_val(ix, iy0) + (1.0 - fy) * ds * p, ix, iy0);
                backproj.set_val(backproj.get_val(ix, iy1) + fy * ds * p, ix, iy1);

                tEnter = t1;
            }
        }
        else
        {
            if (std::abs(dir.y) < 1e-12)
                return;

            const double dt = dVoxelY / std::abs(dir.y);

            while (tEnter < tExit)
            {
                const double t0 = tEnter;
                const double t1 = std::min(t0 + dt, tExit);

                const double tm = 0.5 * (t0 + t1);

                const double x = source.x + tm * dir.x;
                const double y = source.y + tm * dir.y;

                const double px = x / dVoxelX + cx;
                const double py = y / dVoxelY + cy;

                int iy = static_cast<int>(std::floor(py + 0.5));
                iy = std::clamp(iy, 0, ny - 1);

                double pxClamped = std::clamp(
                    px, 0.0, static_cast<double>(nx - 1));

                const int ix0 = static_cast<int>(std::floor(pxClamped));
                const int ix1 = std::min(ix0 + 1, nx - 1);

                const double fx = pxClamped - ix0;

                const double ds = (t1 - t0) * rayLength;

                backproj.set_val(backproj.get_val(ix0, iy) + (1.0 - fx) * ds * p, ix0, iy);
                backproj.set_val(backproj.get_val(ix1, iy) + fx * ds * p, ix1, iy);

                tEnter = t1;
            }
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
            vector2d position_source = config.get_source_position(angle); //fan_beam::get_source_position(angle, config);

            for (std::size_t u = 0; u < config.nDetector; ++u) {
                // Get the position of the voxel in mm and the direction vector between the source and current pixel:
                vector2d position_detector_pixel = config.get_detector_pixel_position(u, angle); //fan_beam::get_detector_pixel_position(u, angle, config);
                vector2d dir = position_detector_pixel - position_source;

                // Get the length of the ray
                const double DSP = std::hypot(dir.x, dir.y);

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
            vector2d position_source = config.get_source_position(angle); //fan_beam::get_source_position(angle, config);

            for (std::size_t u = 0; u < config.nDetector; ++u) {
                // Get the position of the voxel in mm and the direction vector between the source and current pixel:
                vector2d position_detector_pixel = config.get_detector_pixel_position(u, angle); //fan_beam::get_detector_pixel_position(u, angle, config);
                vector2d dir = position_detector_pixel - position_source;

                // Get the length of the ray
                const double DSP = std::hypot(dir.x, dir.y);
                const double p = projections.get_val(u,a);

                fan_beam::accumulate_ray_adjoint(backproj, position_source, dir, p, config);

                // Calculate the contribution of the ray
                // fan_beam::trace_ray(position_source, dir, DSP, config, 
                //     [&](int x, int y, double length) {
                //         backproj.set_val(backproj.get_val(x,y) + length*p, x, y);
                //     }
                // );                
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