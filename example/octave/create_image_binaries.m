clc;
clear all;
close all;

project_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(fullfile(project_dir, 'mex'));

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
    pkg load image
end

if ~exist(fullfile(project_dir, 'img'), 'dir')
    mkdir(fullfile(project_dir, 'img'));
end

%% Shepp-Logan phantom
img = phantom();

fid = fopen(fullfile(project_dir, 'img', 'shepp_logan_gt.bin'), 'wb');
fwrite(fid, size(img, 1), 'double');
fwrite(fid, size(img, 2), 'double');
fwrite(fid, img, 'double');
fclose(fid);
