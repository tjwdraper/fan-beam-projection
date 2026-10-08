clc;
clear all;
close all;

addpath("mex")

pkg load image;

% Load image
img = phantom();

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

% Load mex-function
Ax_mex(img, config);



