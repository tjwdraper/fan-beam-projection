#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>
#include <map>
#include <mex.h>

#include "coord2d.hpp"
#include "mxParser.hpp"
#include "Field.hpp"
#include "fan_beam.hpp"


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

    // Load projections
    opticalflow::Image projections(dim(nDetector, nAngles));
    double* tmp = (double*) mxGetPr(prhs[0]);
    opticalflow::image::load_image(tmp, projections);

    // Allocate memory for adjoint (has image dimensions)
    opticalflow::Image backproj(nVoxel);
    backproj.fill(0.0);

    // Calculate backprojection
    fan_beam::Atb(
        backproj,
        projections,
        DSO, DSD,
        nVoxel, sVoxel, dVoxel,
        nDetector, sDetector, dDetector,
        angles
    );

    // Return projection data to matlab workspace
    mwSize dim_image[2] = {nVoxel.x, nVoxel.y};
    plhs[0] = mxCreateNumericArray(2, dim_image, mxDOUBLE_CLASS, mxREAL);
    tmp = (double*) mxGetPr(plhs[0]);
    opticalflow::image::save_image(tmp, backproj);

}