% 此版本与作者的python代码平行，两者功能上完全相同。
clear all;clc
data.I1 = double(rgb2gray(imread('/Users/Zheng/Documents/pure-version/f/0000.png')));%左 reference
data.I2 = double(rgb2gray(imread('/Users/Zheng/Documents/pure-version/f/0001.png')));
mat = load('/Users/Zheng/Documents/pure-version/f/0000.mat');
data.G1 = mat.G;%外参
data.K = mat.K;%内参
mat = load('/Users/Zheng/Documents/pure-version/f/0001.mat');
data.G2 = mat.G;
data.d0 = ones(size(data.I1))*1.8;

minalpha = 0.015;
dual.warps = 15;
dual.iterations = 20;
dual.scalefactor = 0.75;
dual.levels = 100;
dual.minSize = 48;
dual.stop_level = 0;
dual.check = 10;
dual.lambda = 0.005;
dual.ref = 0;
dual.epsilon = 0.001;
dual.minz = 0.8;
dual.maxz = baseline(data.G1,data.G2)/atan(minalpha/2);

depth = coarse2fine_dual(data, dual);
