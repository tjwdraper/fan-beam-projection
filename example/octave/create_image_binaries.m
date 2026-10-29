clc;
clear all;
close all;

isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
    pkg load image
end

%% Shepp-Logan phantom
img = phantom();

fid = fopen('img/shepp_logan_gt.bin', 'wb');
fwrite(fid, img, 'double');
fclose(fid);

