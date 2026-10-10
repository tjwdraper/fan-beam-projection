clc;
clear all;
close all;

addpath("mex")

pkg load image;

% Load image
x = phantom();
x = x / sqrt(sum(x(:).^2));

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

config.angles = 0:1:100;

for i = 1:100
    y = Ax_mex(x, config);
    z = Atb_mex(y, config);

    x = z ./ sqrt(sum(z(:).^2));

    fprintf("Iter: %d - ||K||: %.5f\n", i, sqrt(sum(x(:).*z(:))));
end



