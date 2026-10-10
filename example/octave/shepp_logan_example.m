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

config.angles = 0:1:100;

% Create projection data
tic;
proj = Ax_mex(img, config);
time = toc;

figure(); imagesc(proj);

% Add TGV parameters to configuration
config.alpha0 = 2.0;
config.alpha1 = 1.0;
config.tau = 0.2;
config.sigma = 0.2;
config.lambda = 0.01;
config.niter = 10;
config.convergence = 1e-6;


img_recon = fan_beam_reconstruction(proj, config);

% Show
figure();
subplot(121); imagesc(img);
subplot(122); imagesc(img_recon);

% Clear functions
clear functions;
