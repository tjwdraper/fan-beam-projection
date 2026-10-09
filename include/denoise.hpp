#ifndef _TGV_DENOISE_2D_HPP_
#define _TGV_DENOISE_2D_HPP_

#include "coord2d.hpp"
#include "Field.hpp"


// Model options
enum class VerboseOption {SILENT, VERBOSE};
enum class ModelOption {TV, TGV};

namespace denoise {
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Norms - used to track convergence:
    ///////////////////////////////////////////////////////////////////////////////////////////////
    double tv_norm(opticalflow::Image& u) {
        const dim dimin = u.get_dimensions();

        double norm(0.0);
        for (std::size_t j = 0; j < dimin.y; ++j) {
            for (std::size_t i = 0; i < dimin.x; ++i) {
                double uxv = opticalflow::gradients::dx(u, i, j);
                double uyv = opticalflow::gradients::dy(u, i, j);

                norm += std::sqrt(uxv*uxv + uyv*uyv);
            }
        }
        return norm;
    }

    double tv_residual_norm(const opticalflow::Image& u, const opticalflow::Image& vx, const opticalflow::Image& vy) {
        const dim dimin = u.get_dimensions();

        double norm(0.0);
        for (std::size_t j = 0; j < dimin.y; ++j) {
            for (std::size_t i = 0; i < dimin.x; ++i) {
                double uxv = opticalflow::gradients::dx(u, i, j) - vx.get_val(i,j);
                double uyv = opticalflow::gradients::dy(u, i, j) - vy.get_val(i,j);

                norm += std::sqrt(uxv*uxv + uyv*uyv);
            }
        }
        return norm;
    }

    double l1_E_norm(const opticalflow::Image& vx, const opticalflow::Image& vy) {
        const dim dimin = vx.get_dimensions();

        double norm(0.0);
        for (std::size_t j = 0; j < dimin.y; ++j) {
            for (std::size_t i = 0; i < dimin.x; ++i) {
                const double Evbarxx = opticalflow::gradients::dx(vx, i, j);
                const double Evbaryy = opticalflow::gradients::dy(vy, i, j);
                const double Evbarxy = 0.5 * (opticalflow::gradients::dx(vy,i,j) + 
                                              opticalflow::gradients::dy(vx,i,j));

                norm += std::sqrt(Evbarxx*Evbarxx + Evbaryy*Evbaryy + 2*Evbarxy*Evbarxy);
            }
        }
        return norm;

    }

    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Update operators
    ///////////////////////////////////////////////////////////////////////////////////////////////
    void update_p(opticalflow::Image& px, opticalflow::Image& py, const opticalflow::Image& ubar, double sigma) {        
        const dim dimin = px.get_dimensions();

        for (std::size_t j = 0; j < dimin.y; ++j) {
            for (std::size_t i = 0; i < dimin.x; ++i) {
                px.set_val(px.get_val(i,j) + sigma * (opticalflow::gradients::dx(ubar, i, j)),i,j);
                py.set_val(py.get_val(i,j) + sigma * (opticalflow::gradients::dy(ubar, i, j)),i,j);
            }
        }
    }
    
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

    void update_r(opticalflow::Image& r, const opticalflow::Image& ubar, const opticalflow::Image& p, double sigma) {
        
    }

    void update_u(opticalflow::Image& u,
                  const opticalflow::Image& px, const opticalflow::Image& py,
                  double tau) {
        const dim dimin = u.get_dimensions();        

        for (std::size_t j = 0; j < dimin.y; ++j) {
            for (std::size_t i = 0; i < dimin.x; ++i) {
                const double divp = opticalflow::gradients::adjoint_dx(px,i,j) + opticalflow::gradients::adjoint_dy(py,i,j);
                u.set_val(u.get_val(i,j) + tau*divp,i,j);
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

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Projection and proximal operators
    ///////////////////////////////////////////////////////////////////////////////////////////////
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

    void prox(opticalflow::Image& u, const opticalflow::Image& f, double tau, double lambda) {
        const dim dimin = u.get_dimensions();

        for (std::size_t j = 0; j < dimin.y; ++j) {
            for (std::size_t i = 0; i < dimin.x; ++i) {
                u.set_val( (lambda*u.get_val(i,j) + tau * f.get_val(i,j)) / (lambda + tau), i, j);
            }
        }
    }

    // The primal-dual TGV-denoising algorithm
    opticalflow::Image tgv_reconstruct(const opticalflow::Image& p, 
        double tau, double sigma, double lambda, double alpha0, double alpha1, int niter, double convergence, 
        double DSO, double DSD, dim nVoxel, vector2d sVoxel, vector2d dVoxel, int nDetector, double sDetector, double dDetector, const std::vector<double>& angles) {
        
        const int nAngles = angles.size();

        // Initialize Primal-Dual variables
        opticalflow::Image u(nVoxel); u.fill(0.0);
        opticalflow::Image ubar(nVoxel); ubar.fill(0.0);

        opticalflow::Image vx(nVoxel); vx.fill(0.0);
        opticalflow::Image vy(nVoxel); vy.fill(0.0);

        opticalflow::Image vbarx(nVoxel); vbarx.fill(0.0);
        opticalflow::Image vbary(nVoxel); vbary.fill(0.0);

        opticalflow::Image px(nVoxel); px.fill(0.0);
        opticalflow::Image py(nVoxel); py.fill(0.0);

        opticalflow::Image qxx(nVoxel); qxx.fill(0.0);
        opticalflow::Image qyy(nVoxel); qyy.fill(0.0);
        opticalflow::Image qxy(nVoxel); qxy.fill(0.0);

        opticalflow::Image r(dim(nDetector, nAngles)); r.fill(0.0);

        // Tracking variables
        opticalflow::Image uold(nVoxel); uold.fill(0.0);
        opticalflow::Image vxold(nVoxel); vxold.fill(0.0);
        opticalflow::Image vyold(nVoxel); vyold.fill(0.0);    

        // Primal-dual iterations
        for (int iter = 0; iter < niter; ++iter) {
            // Update dual variables
            denoise::update_p(px,py,ubar,vbarx,vbary,sigma);
            denoise::update_q(qxx, qyy, qxy, vbarx, vbary, sigma);
            denoise::update_r(r, ubar, p, sigma);

            denoise::proj_p(px,py,alpha1);
            denoise::proj_q(qxx,qyy,qxy,alpha0);
            denoise::prox_r(r, sigma, lambda);

            // Track u
            uold = u;

            // Proximal operator
            denoise::update_u(u, px, py, r, tau);

            // Update ubar
            ubar = 2*u - uold;

            // Track v
            vxold = vx;
            vyold = vy;

            // Update v
            denoise::update_v(vx, vy, px, py, qxx, qyy, qxy, tau);

            // Update vbar
            vbarx = 2*vx - vxold;
            vbary = 2*vy - vyold;

            // Check for convergence
            double relchange = opticalflow::image::norm(u-uold)/opticalflow::image::norm(u);
            if (relchange < convergence)
                break;

            // Update some norms:
            // if (iter % 50 == 0) {
            //     std::cout << "iter " << iter 
            //         << "\t|u_new-u|/|u_new| = " << relchange 
            //         << "\t1/(2*lambda)*|u-f|^2 = " << opticalflow::image::normsq(u-f) / (2*lambda)
            //         << "\talpha_1 * |Du - v|_1 = " << denoise::tv_residual_norm(u,vx,vy) * alpha1
            //         << "\talpha_0 * |Ev|_1 = " << denoise::l1_E_norm(vx,vy) * alpha0
            //         << std::endl;
            // }
        }

        return u;
    }


    opticalflow::Image tgv_denoise(const opticalflow::Image& f, double tau, double sigma, double lambda, double alpha0, double alpha1, int niter, double convergence, VerboseOption verbose) {
        if (tau <= 0.0)
            throw std::runtime_error("Tau has to be a positive scalar.");

        if (sigma <= 0) 
            throw std::runtime_error("Sigma has to be a positive scalar.");

        // Get image dimensions
        const dim dimin = f.get_dimensions();

        // Initialize Primal-Dual variables
        opticalflow::Image u(f);
        opticalflow::Image ubar(f);

        opticalflow::Image vx(dimin); vx.fill(0.0);
        opticalflow::Image vy(dimin); vy.fill(0.0);

        opticalflow::Image vbarx(dimin); vbarx.fill(0.0);
        opticalflow::Image vbary(dimin); vbary.fill(0.0);

        opticalflow::Image px(dimin); px.fill(0.0);
        opticalflow::Image py(dimin); py.fill(0.0);

        opticalflow::Image qxx(dimin); qxx.fill(0.0);
        opticalflow::Image qyy(dimin); qyy.fill(0.0);
        opticalflow::Image qxy(dimin); qxy.fill(0.0);

        // Tracking variables
        opticalflow::Image uold(dimin); uold.fill(0.0);
        opticalflow::Image vxold(dimin); vxold.fill(0.0);
        opticalflow::Image vyold(dimin); vyold.fill(0.0);

        // Primal-dual iterations
        for (int iter = 0; iter < niter; ++iter) {
            // Update dual variables
            denoise::update_p(px,py,ubar,vbarx,vbary,sigma);
            denoise::update_q(qxx, qyy, qxy, vbarx, vbary, sigma);

            denoise::proj_p(px,py,alpha1);
            denoise::proj_q(qxx,qyy,qxy,alpha0);

            // Track u
            uold = u;

            // Proximal operator
            denoise::update_u(u, px, py, tau);
            denoise::prox(u, f, tau, lambda);

            // Update ubar
            ubar = 2*u - uold;

            // Track v
            vxold = vx;
            vyold = vy;

            // Update v
            denoise::update_v(vx, vy, px, py, qxx, qyy, qxy, tau);

            // Update vbar
            vbarx = 2*vx - vxold;
            vbary = 2*vy - vyold;

            // Check for convergence
            double relchange = opticalflow::image::norm(u-uold)/opticalflow::image::norm(u);
            if (relchange < convergence)
                break;

            // Update some norms:
            if (iter % 50 == 0 && verbose == VerboseOption::VERBOSE) {
                std::cout << "iter " << iter 
                    << "\t|u_new-u|/|u_new| = " << relchange 
                    << "\t1/(2*lambda)*|u-f|^2 = " << opticalflow::image::normsq(u-f) / (2*lambda)
                    << "\talpha_1 * |Du - v|_1 = " << denoise::tv_residual_norm(u,vx,vy) * alpha1
                    << "\talpha_0 * |Ev|_1 = " << denoise::l1_E_norm(vx,vy) * alpha0
                    << std::endl;
            }
        }

        return u;
    }

    
}

#endif