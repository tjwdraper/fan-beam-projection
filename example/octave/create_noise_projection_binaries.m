clc;
clear all;
close all;

dir_name = fileparts(mfilename('fullpath'));
addpath([dir_name, '\..\..\mex']);

%% Select I0
I0 = 1e4;

%% Check if script is launched from Matlab or GNU Octave
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
    pkg load image;
    pkg load statistics;
end

%% Read fan beam geometry and
json_data = jsondecode(fileread('./example/config_shepp_logan.json'));

geometry = json_data.geometry;
geometry.nVoxel = int32(geometry.nVoxel);
if (geometry.nAngles == 1 && geometry.minAngle == geometry.maxAngle)
    geometry.angles = geometry.minAngle;
elseif (geometry.nAngles > 1 && geometry.minAngle < geometry.maxAngle)
    minAngle = geometry.minAngle;
    maxAngle = geometry.maxAngle;
    dAngle = (maxAngle - minAngle) / (geometry.nAngles - 1);
    geometry.angles = minAngle:dAngle:maxAngle;
else
    error("invalid angle configuration given in json file.");
end

%% Read image binary. If not in the image folder, create the image binary
if ~exist([dir_name, '\..\..\img'], 'dir')
  mkdir('img');
end

if ~exist([dir_name, '\..\..\img\shepp_logan_gt.bin'], 'file')
  create_image_binaries;
end

fid = fopen([dir_name, '\..\..\img\shepp_logan_gt.bin'], 'rb');
nx = fread(fid, 1, 'double');
ny = fread(fid, 1, 'double');
if (nx ~= geometry.nVoxel(1) || ...
    ny ~= geometry.nVoxel(2))
    fclose(fid);
    error('image dimensions read from binary do not match with json configuration.');
end
img = fread(fid, geometry.nVoxel, 'double');
fclose(fid);

%% Apply forward projection operator
projections = Ax_mex(img, geometry);

fid = fopen([dir_name, '\..\..\img\projections_gt.bin'], 'wb');
fwrite(fid, size(projections, 1), 'double');
fwrite(fid, size(projections, 2), 'double');
fwrite(fid, projections, 'double');
fclose(fid);

%% Add Poisson noise

% Convert projection to photon count, add Poisson noise and save result
I = I0 * exp(-projections / max(projections(:)));
Inoise = poissrnd(I);
projections_noise =  -max(projections(:)) * log(Inoise / I0);

fid = fopen([dir_name, '\..\..\img\projections_noise.bin'], 'wb');
fwrite(fid, size(projections_noise, 1), 'double');
fwrite(fid, size(projections_noise, 2), 'double');
fwrite(fid, projections_noise, 'double');
fclose(fid);

fprintf("MSE: %.3f\n", sum(projections(:) - projections_noise(:)).^2)
figure();
subplot(121); imagesc(projections);
subplot(122); imagesc(projections_noise);


