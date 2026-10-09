#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <map>
#include <mex.h>
#include <chrono>

#include "coord2d.hpp"
#include "mxParser.hpp"
#include "Field.hpp"
#include "denoise.hpp"

///////////////////////////////////////////////////////////////////////////////////////////////////
// Main function
///////////////////////////////////////////////////////////////////////////////////////////////////
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if ((nlhs != 1) || (nrhs != 2))
        mexErrMsgTxt("Error: invalid number of input and output variables gives.\n");

    // Read geometry
    double DSO; mxParser::parse_DSO(DSO, prhs[1]);
    double DSD; mxParser::parse_DSD(DSD, prhs[1]);

    dim nVoxel; mxParser::parse_nVoxel(nVoxel, prhs[1]);
    vector2d sVoxel; mxParser::parse_sVoxel(sVoxel, prhs[1]);
    vector2d dVoxel; mxParser::parse_dVoxel(dVoxel, prhs[1]);

    int nDetector; mxParser::parse_nDetector(nDetector, prhs[1]);
    double dDetector; mxParser::parse_dDetector(dDetector, prhs[1]);
    double sDetector; mxParser::parse_sDetector(sDetector, prhs[1]);

    std::vector<double> angles;
    mxParser::parse_angles(angles, prhs[1]);

    int nAngles = angles.size();

    // Read the TGV parameters
    double alpha0; mxParser::parse_alpha0(alpha0, prhs[1]);
    double alpha1; mxParser::parse_alpha1(alpha1, prhs[1]);
    double tau; mxParser::parse_tau(tau, prhs[1]);
    double sigma; mxParser::parse_sigma(sigma, prhs[1]);
    double lambda; mxParser::parse_lambda(lambda, prhs[1]);
    int niter; mxParser::parse_niter(niter, prhs[1]);
    double convergence; mxParser::parse_convergence(convergence, prhs[1]);

    // Load projections
    opticalflow::Image projections(dim(nDetector, nAngles));
    double* tmp = (double*) mxGetPr(prhs[0]);
    opticalflow::image::load_image(tmp, projections);

    // Allocate memory for reconstructed image
    opticalflow::Image image(nVoxel);
    image.fill(0.0);

    // Reconstruct image
    const auto start = std::chrono::steady_clock::now();
    image = denoise::tgv_reconstruct(
        projections, 
        tau, sigma, lambda, alpha0, alpha1, niter, convergence,
        DSO, DSD, nVoxel, sVoxel, dVoxel, nDetector, sDetector, dDetector, angles
    );
    const auto end = std::chrono::steady_clock::now();

    const std::chrono::duration<double> elapsed = end-start;
    std::cout << "Runtime (s): " << elapsed.count() << "\n";

    // Return image to matlab workspace
    mwSize dim_image[2] = {nVoxel.x, nVoxel.y};
    plhs[0] = mxCreateNumericArray(2, dim_image, mxDOUBLE_CLASS, mxREAL);
    tmp = (double*) mxGetPr(plhs[0]);
    opticalflow::image::save_image(tmp, image);
}