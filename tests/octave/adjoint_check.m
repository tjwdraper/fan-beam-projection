clc;
clear all;
close all;

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
    pkg load image;
end

addpath("mex")

% Set fan-beam geometry
config = struct();
config.DSO = 1000;
config.DSD = 1536;

config.nVoxel = int32([256, 256]);
config.sVoxel = [256, 256];
config.dVoxel = [1, 1];

config.nDetector = int32(512);
config.dDetector = 0.8;
config.sDetector = 409.6;

config.angles = 0:10:100;

% Create some artificial values - image and projections
p = randn(config.nDetector, length(config.angles));
m = randn(config.nVoxel(1), config.nVoxel(2));

Am = Ax_mex(m, config);
Atp = Atb_mex(p, config);

ptAm = sum(p(:).*Am(:));
Atpm = sum(Atp(:).*m(:));

% A is self-adjoint if <A^Tp, m> = <p, Am>
fprintf("Self adjointness: <A^Tp,m> - <p,Am> = %.12f\n", ptAm - Atpm);

clear functions;
