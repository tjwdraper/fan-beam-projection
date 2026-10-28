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

    // Read fan-beam configuration
    fan_beam_configuration config = mxParser::parse_fan_beam_configuration(prhs[1]);

    // Load image
    opticalflow::Image image(config.nVoxel);
    double* tmp = (double*) mxGetPr(prhs[0]);
    opticalflow::image::load_image(tmp, image);

    // Allocate memory for projection data
    opticalflow::Image projections(dim(config.nDetector, config.nAngles));
    projections.fill(0.0);

    // Calculate projections
    fan_beam::Ax(
        projections,
        image,
        config
    );

    // Return projection data to matlab workspace
    mwSize dim_projections[2] = {config.nDetector, config.nAngles};
    plhs[0] = mxCreateNumericArray(2, dim_projections, mxDOUBLE_CLASS, mxREAL);
    tmp = (double*) mxGetPr(plhs[0]);
    opticalflow::image::save_image(tmp, projections);

}