Copyright (c) ICG. All rights reserved.

Institute for Computer Graphics and Vision
Graz University of Technology / Austria

CONTACT:
Gottfried Graber: graber@icg.tugraz.at

=============================================================================

This is a python implementation of the minimal surface regularization of 
perspective depthmaps. It has been written and tested under Linux using
python2.7, although it should be straightforward to use it in Windows and/or 
python3.x

Usage: Load flowdepth.py into your python IDE and run it, or start it directly
from the terminal via
$ python flowdepth.py

Sample input data is provided, feel free to try it with your own data! The
framework expects a dictionary which contains the following:
*) a key 'images', value: a list of 2 numpy-images, first one is the reference image
*) a key 'G', value: list with the camera transformation matrices for the images
*) a key 'K', value: intrinsic calibration matrix
See line 90 how the dictionary is constructed.

=============================================================================

If you use this code, please cite the following publication:

@inproceedings{graber_cvpr2015,
author = {Gottfried Graber and Jonathan Balzer and Stefano Soatto and Thomas Pock},
title = {{Efficient Minimal-Surface Regularization of Perspective Depth Maps in Variational Stereo}},
booktitle = {IEEE Conference on Computer Vision and Pattern Recognition},
year = {2015}
}

For further questions, feel free to contact via: graber@icg.tugraz.at