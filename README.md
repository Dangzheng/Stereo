# Stereo
The repository contains code about stereo and some practise. Dataset generate code made by myself.

## matlab-test-dataset

This folder contains matlab code for generating synthetic plane, slanted plane as foreground, background plane, intrinsic parametes, external parameters and groundtruth.

## matlab-l2norm-minimun surface

This part of code using l2-norm and gradient descent，the regular item of the perspective model realizing by minimizing the surface area was used in this paper. While I utilized in this place is just l2-norm,not l2，1-norm in the paper，for the simply implementation of l2-norm, which can be solved by Gradient Descent Method，instead of primal-dual. This can help test for correctness of other place in code. The reason why I choose to rewrite to the version of matlab is that it is so convenient to debug.

**Efficient Minimal-Surface Regularizationof Perspective Depth Maps in Variational Stereo**

## matlab-minimum surface

This code is corresponding to **Efficient Minimal-Surface Regularizationof Perspective Depth Maps in Variational Stereo**

I just rewrite it by matlab.
## python-l0smooth
This code is corresponding to **Image Smoothing via L0 Gradient Minimization**. I use python to rewrite it.

## myself-dataset-test-minimal-surface
This fold loads the pictures producting by my own，testing on code minimal-surface，using the author's python code，editting in pycharm。

The final reason for the problem is normalization problem，the operating range in author's code is [0.8, 2.48]，as for the reason why it can only be calculated in such small depth range, it needed for further study.

## middlebury_minimal_surface
This fold loads the python code either，which is the non GPU version published by author. The test dataset is Middlebury Stereo 2014. The unit using in this dataset is the ddcimal inch. I should try to convert it. At the same time it is confused that why it only can be calculated in small range in author's code.
