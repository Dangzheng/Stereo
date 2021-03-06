clear;clc
data.I1 = double(rgb2gray(imread('/Users/Zheng/Documents/MATLAB/venus/im2.ppm')));%左 reference
data.I2 = double(rgb2gray(imread('/Users/Zheng/Documents/MATLAB/venus/im6.ppm')));
[r,c] = size(data.I1);
baseline = 2;
data.G = [1 0 0 -baseline;0 1 0 0;0 0 1 0];%外参
data.K = [c 0 c/2;0 c r/2;0 0 1];%内参
%data.d0 = ones(size(data.I1))*1.8;
disp = double(imread('/Users/Zheng/Documents/MATLAB/venus/disp2.pgm') );
%data.d0 = disp2depth(disp, c, baseline);
data.d0 = disp;
minalpha = 0.015;
params.alpha = 1;%数据项系数
params.beta = 0.1;%link项系数
params.betamax = 1e5;
params.kappa = 2;
params.lambda = 2;%l0平滑项系数
params.minz = 95;
params.maxz = baseline/atan(minalpha/2);
params.scalefactor = 0.95;
params.level_max = 100;
params.stop_level = 1;
params.minSize = 48;
params.itermax = 50;
params.epsilon = 0.001;

dual.warps = 15;
dual.iterations = 20;
dual.scalefactor = params.scalefactor;
dual.levels = 100;
dual.minSize = params.minSize;
dual.stop_level = 0;
dual.check = 10;
dual.lambda = 0.5;
dual.ref = 0;
dual.epsilon = 0.001;
dual.minz = 0.8;
dual.maxz = params.maxz;
%depth = coarse2fine(data, params);
depth = coarse2fine_dual(data, params, dual);
