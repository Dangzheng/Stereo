# -*- coding: utf-8 -*-
"""
Created on Wed Feb  4 16:34:56 2015

@author: gottfried
"""

import numpy as np
from scipy import ndimage

from flowtools import to_z
from flowtools import backproject
import time

    
def warp_zeta(Iref, I, Grel, K, zeta):
    Xn = backproject(zeta.shape, K)
    #backproject 是利用相机的内参矩阵做变换
    X = np.vstack((Xn, 1/to_z(zeta.flatten())))    # homogeneous coordinate = inverse depth
    #相当于就是在最后一行上面都加上了每个点的深度的倒数
    #print(Grel.shape)3x4
    if Grel.shape[0] < 4:
        Grel = np.vstack((Grel, np.array([0,0,0,1])))
    #进行补行之后Grel为4x4    
    dX = np.dot(Grel[0:3,0:3], Xn)*1/to_z(zeta.flatten())   # derivative of transformation G*X
    
    X2 = np.dot(Grel, X)             # transform point
    X2 = np.dot(Grel, X)             # transform point
    x2 = np.dot(K, X2[0:3,:])        # project to image plane
    
    x2[0,:] /= x2[2,:]             # dehomogenize
    x2[1,:] /= x2[2,:]
    Iw = ndimage.map_coordinates(I, np.vstack(np.flipud(x2[0:2,:])), order=1, cval=np.nan)
    #flipud是实现矩阵的上下翻转
    Iw = np.reshape(Iw, I.shape)      # warped image
    gIwy,gIwx = np.gradient(Iw)
    #这一步求x,y方向上面的梯度，对应的是文章中(38)式
    fx = K[0,0]; fy = K[1,1]
    for i in range(0,3):           # dehomogenize
        X2[i,:] /= X2[3,:]
    z2 = X2[2,:]**2
    dT = np.zeros((2,X2.shape[1]))
    dT[0,:] = fx/X2[2,:]*dX[0,:] - fx*X2[0,:]/z2*dX[2,:]    # derivative of projected point x2 = pi(G*X)
    dT[1,:] = fy/X2[2,:]*dX[1,:] - fy*X2[1,:]/z2*dX[2,:]
    
    # full derivative I(T(x,z)), T(x,z)=pi(G*X)
    Ig = np.reshape(gIwx.flatten()*dT[0,:] + gIwy.flatten()*dT[1,:], zeta.shape)
    It = Iw - Iref       # 'time' derivative'
    
    It[np.where(np.isnan(Iw))] = 0;
    Ig[np.where( (np.isnan(gIwx).astype('uint8') + np.isnan(gIwy).astype('uint8'))>0 ) ] = 0
    
    return Iw, It, Ig
