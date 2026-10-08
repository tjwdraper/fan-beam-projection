#ifndef _MX_PARSER_HPP_
#define _MX_PARSER_HPP_

#include <string>
#include <mex.h>
#include <type_traits>
#include <optional>
#include <algorithm>
#include <cstddef>
#include <map>

#include "coord2d.hpp"
#include "ConfigurationOptions.hpp"

namespace mxParser {
    inline VerboseOption verbose = VerboseOption::VERBOSE;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Helper functions
///////////////////////////////////////////////////////////////////////////////////////////////////
template <typename T>
T parse_scalar(const mxArray* config, const char* fieldname, const std::optional<T> default_value = std::nullopt) {
    const mxArray* field = mxGetField(config, 0, fieldname);
    if (!field && !default_value.has_value()) {
        mexErrMsgIdAndTxt("mxParser:InputError", "field %s not found.", fieldname);
    }
    else if (!field && default_value.has_value()) {
        if (mxParser::verbose == VerboseOption::VERBOSE) {
            mexWarnMsgIdAndTxt("mxParser:InputWarning", "field %s not found. set to default value.", fieldname);
        }
        return default_value.value();
    }

    if (!mxIsScalar(field)) {
        mexErrMsgIdAndTxt("mxParser:TypeError", "field %s must be a scalar.", fieldname);
    }

    if (mxIsComplex(field)) {
        mexErrMsgIdAndTxt("mxParser:TypeError", "field %s must contain real values.", fieldname);
    }
    return static_cast<T>(mxGetScalar(field));
}

template <typename T>
T parse_string(const mxArray* config, const char* fieldname, const std::map<std::string, T>& mapper, const std::optional<T> default_value = std::nullopt) {
    const mxArray* field = mxGetField(config, 0, fieldname);
    if (!field && !default_value.has_value()) {
        mexErrMsgIdAndTxt("mxParser:InputError", "field %s not found.", fieldname);
    }
    else if (!field && default_value.has_value()) {
        if (mxParser::verbose == VerboseOption::VERBOSE) {
            mexWarnMsgIdAndTxt("mxParser:InputWarning", "field %s not found. set to default value.", fieldname);        
        }
        return default_value.value();
    }

    if (!mxIsChar(field)) {
        mexErrMsgIdAndTxt("mxParser:TypeError", "field %s is not a valid character array (did you try \" instead of \'?).", fieldname);
    }

    char* field_c = mxArrayToString(field);
    std::string key(field_c);
    mxFree(field_c);

    auto it = mapper.find(key);
    if (it == mapper.end()) {
        mexErrMsgIdAndTxt("mxParser:InputError", "invalid field option given for field name %s.", fieldname);
    }

    return it->second;
}

inline dim parse_dim(const mxArray* config, const char* fieldname, const std::optional<dim> default_value = std::nullopt) {
    const mxArray* field = mxGetField(config, 0, fieldname);
    if (!field && !default_value.has_value()) {
        mexErrMsgIdAndTxt("mxParser:InputError", "field %s not found.", fieldname);
    }
    else if (!field && default_value.has_value()) {
        if (mxParser::verbose == VerboseOption::VERBOSE) {
            mexWarnMsgIdAndTxt("mxParser:InputWarning", "field %s not found. set to default value.", fieldname);        
        }
        return default_value.value();
    }

    if (!mxIsInt32(field)) {
        mexErrMsgIdAndTxt("mxParser:TypeError", "field %s is not a valid int32 array.", fieldname);
    }

    if (mxIsComplex(field)) {
        mexErrMsgIdAndTxt("mxParser:TypeError", "field %s must contain real values.", fieldname);
    }

    if (mxGetNumberOfElements(field) != 2) {
        mexErrMsgIdAndTxt("mxParser:TypeError", "field %s is not of length two (2).", fieldname);
    }

    int32_T* dim_arr = static_cast<int32_T*>(mxGetData(field));
    return dim(dim_arr[0], dim_arr[1]);
}

template <typename T>
T* parse_array_ptr(const mxArray* config, const char* fieldname) {
    const mxArray* field = mxGetField(config, 0, fieldname);
    if (!field) {
         mexErrMsgIdAndTxt("mxParser:InputError", "field %s not found.", fieldname);        
    }
    if (!mxIsNumeric(field)) {
        mexErrMsgIdAndTxt("mxParser:TypeError", "field %s does not contain numeric data.", fieldname);
    }
    if (mxIsComplex(field)) {
        mexErrMsgIdAndTxt("mxParser:TypeError", "field %s must contain real values.", fieldname);
    }

    if constexpr (std::is_same_v<T, int32_T>) {
        if (mxGetClassID(field) != mxINT32_CLASS) {
            mexErrMsgIdAndTxt("mxParser:TypeError", "the array contained in field %s must be type int32.", fieldname);
        }
        return static_cast<int32_T*>(mxGetData(field));
    }
    else if constexpr (std::is_same_v<T, double>) {
        if (mxGetClassID(field) != mxDOUBLE_CLASS) {
            mexErrMsgIdAndTxt("mxParser:TypeError", "the array contained in field %s must be type double.", fieldname);
        }
        return static_cast<double*>(mxGetData(field));
    }
    else {
        mexErrMsgIdAndTxt("mxParser:TypeError", "the array contained in field %s must be either of type int32 or double.", fieldname);
    }
}

inline std::size_t parse_array_size(const mxArray* config, const char* fieldname) {
    const mxArray* field = mxGetField(config, 0, fieldname);
    if (!field) {
         mexErrMsgIdAndTxt("mxParser:InputError", "field %s not found.", fieldname);        
    }
    if (!mxIsNumeric(field)) {
        mexErrMsgIdAndTxt("mxParser:TypeError", "field %s does not contain numeric data.", fieldname);
    }
    return static_cast<std::size_t>(mxGetNumberOfElements(field));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// mxParser namespace
///////////////////////////////////////////////////////////////////////////////////////////////////
namespace mxParser {
    // inline VerboseOption verbose = VerboseOption::VERBOSE;

    // Generic method: check if field exists in Matlab/GNU Octave structure
    inline bool exists(const mxArray* config, const char* fieldname) {
        const mxArray* field = mxGetField(config, 0, fieldname);
        return field != nullptr;
    }

    // Parse verbose options
    inline void parse_verbose_option(const mxArray* config) {
        mxParser::verbose = parse_string<VerboseOption>(config, "verbose_option", mapper_verbose_option, VerboseOption::VERBOSE);
    }

    // Parse optical flow method
    inline void parse_opticalflow_option(ModelOption& option, const mxArray* config) {
        option = parse_string<ModelOption>(config, "optical_flow_option", mapper_model_option, ModelOption::HORN_SCHUNCK);
    }

    // Parse registration parameters
    inline void parse_size_image(dim& dimin, const mxArray* config) {
        dimin = parse_dim(config, "size_image");
    }
    inline void parse_nscales(std::size_t& nscales, const mxArray* config) {
        nscales = parse_array_size(config, "niter") - 1;
    }
    inline void parse_niter(std::size_t* niter, const mxArray* config) {
        std::size_t nscales;
        mxParser::parse_nscales(nscales, config);

        int32_T* niter_arr = parse_array_ptr<int32_T>(config, "niter");
        for (std::size_t i = 0; i < nscales + 1; ++i) {
            if (niter_arr[i] < 0)
                mexErrMsgIdAndTxt("mxParser:ValueError", "field %s contains negative values.", "niter");
            niter[i] = static_cast<std::size_t>(niter_arr[i]);
        }
    }
    inline void parse_alpha(double& alpha, const mxArray* config) {
        alpha = parse_scalar<double>(config, "alpha", std::optional<double>(0.4));
        if (alpha < 0.0)
                mexErrMsgIdAndTxt("mxParser:ValueError", "alpha has to be positive");
    }
    inline void parse_beta(double& beta, const mxArray* config) {
        beta = parse_scalar<double>(config, "beta", std::optional<double>(0.2));
        if (beta < 0.0)
                mexErrMsgIdAndTxt("mxParser:ValueError", "beta has to be positive.");
    }
    inline void parse_convergence_threshold(double& eps, const mxArray* config) {
        eps = parse_scalar<double>(config, "eps", std::optional<double>(1e-3));
    }
    inline void parse_nrefine(std::size_t& nrefine, const mxArray* config) {
        nrefine = parse_scalar<std::size_t>(config, "nrefine", std::optional<std::size_t>(0));
    }

    inline void parse_resampling_factor(double& resampling_factor, const mxArray* config) {
        resampling_factor = parse_scalar<double>(config, "resampling_factor", std::optional<double>(0.5));
        if (resampling_factor <= 0 || resampling_factor >= 1)
                mexErrMsgIdAndTxt("mxParser:ValueError", "field resampling_factor has to contain a value between 0 and 1.");
    }
}

#endif