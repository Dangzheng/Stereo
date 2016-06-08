clear;clc
data.I1 = double(rgb2gray(imread('/Users/Zheng/Documents/MATLAB/venus/im2.ppm')));%�� reference
data.I2 = double(rgb2gray(imread('/Users/Zheng/Documents/MATLAB/venus/im6.ppm')));
[r,c] = size(data.I1);
baseline = 2;
data.G = [1 0 0 -baseline;0 1 0 0;0 0 1 0];%���
data.K = [c 0 c/2;0 c r/2;0 0 1];%�ڲ�

input = 'gauss';
if input == 'gauss'
   disp = double( imread('/Users/Zheng/Documents/MATLAB/venus/disp2.pgm') );
   init = disp2depth(disp);
   data.d0 = randn(size(init))*0 + init;
else
   data.d0 = ones(size(data.I1))*50;
end

minalpha = 0.015;
params.alpha = 1;%������ϵ��
params.beta = 2;%link��ϵ��
params.betamax = 1e5;
params.kappa = 2;
params.lambda = 0.0005;%l0ƽ����ϵ��
params.minz = 95;
params.maxz = baseline/atan(minalpha/2);
params.scalefactor = 0.95;
params.level_max = 100;
params.stop_level = 1;
params.minSize = 48;
params.itermax = 50;
params.epsilon = 0.001;


depth = coarse2fine(data, params);