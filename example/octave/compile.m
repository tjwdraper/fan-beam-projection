clc;

#mkoctfile --mex -O3 -Iinclude/ -o ./mex/Ax_mex.mex ./src/Ax_mex.cpp
#mkoctfile --mex -O3 -Iinclude/ -o ./mex/Atb_mex.mex ./src/Atb_mex.cpp
mkoctfile --mex -O3 -Iinclude/ -o ./mex/fan_beam_reconstruction.mex ./src/fan_beam_reconstruct.cpp
