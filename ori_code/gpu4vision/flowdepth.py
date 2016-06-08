# -*- coding: utf-8 -*-
"""
Created on Tue Feb  3 16:29:58 2015

@author: gottfried
"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt
from PyQt4 import QtGui
import scipy.io as sio

import flowtools as ft
import warping as warp
import solver as pd
import cv2
import time

def run_ctf(params, d0, data, imageWidget=None):
    #给进来的参数，depth map的初始值，数据（图片，G，K） 
    levels = (np.floor( (np.log(params['minSize'])-np.log(np.minimum(*d0.shape))) / np.log(params['scalefactor']) ) + 1).astype('int')
    levels = np.minimum(levels, params['levels'])
    warps = params['warps']#为什么会直接就设置了一个warps
    iterations = params['iterations']
    
    L = params['Lambda']          # compute lambda sequence for the levels
    Lambda = np.zeros(levels)
    for idx,l in enumerate(Lambda):
        Lambda[idx] = L*(1.0/params['scalefactor'])**idx#此处的data term前面的系数lambda 就是我的alpha
    #Lambda 系数还是随着迭代的level不同而在不断的变化
 
    params['Lambda'] = Lambda
    dim_dual = 3        # dimension of dual variable
    
    for level in range(levels-1,params['stop_level']-1,-1):
        level_sz = np.round(np.array(d0.shape) * params['scalefactor']**level)
        factor = np.hstack(( level_sz / np.array(d0.shape), 1 ))
        K = data['K'] * (np.ones((3,3))*factor).T; K[2,2] = 1    # scale K matrix according to image size
                                                               # individual scaling for x & y
        print '--- level %d, size %.1f %.1f' % (level, level_sz[1], level_sz[0])
        
        if level == levels-1:
            d = ft.imresize(d0, level_sz)      # initialization at coarsest level
            p = np.zeros(d.shape + (dim_dual,))
        else:
            print 'prolongate'
            d = ft.imresize(d, level_sz)          # prolongate to finer level
            ptmp = p.copy()
            p = np.zeros(d.shape+(dim_dual,))
            for i in range(0, dim_dual):
                p[:,:,i] = ft.imresize(ptmp[:,:,i], level_sz)
                
        Xn = ft.backproject(d.shape, K)#相当于将点都利用相机内参，变换到了perspective parameterization
        L_normal = ft.make_linearOperator(d.shape, Xn, K)
        print(L_normal);time.sleep(9999)
        #对应的是文章中（21）的公式L而且表示成了稀疏的形式。
                
        img_scaled = []                    # scale all images
        for img in data['images']:
            img_scaled.append(ft.imresize(img, level_sz))

        
        Iref = img_scaled[0]#此处用的是缩小size之后的图像
        Gref = data['G'][0]#此处的G为文章中提及的两个pi中间夹的g的取值。4x4
        
            
        for k in range(1,warps+1):
            print 'warp', k, 'maxZ=', params['maxz']
            
            warpdata = []                       # warp all images
            for img,g in zip(img_scaled[1:], data['G'][1:]):
                Grel = ft.relative_transformation(Gref, g)
                #此处的Grel是左右两个相机之间的相对关系，就是左相机转换到右相机坐标系的转换矩阵
                print(d)
                warpdata.append(warp.warp_zeta(Iref, img, Grel, K, ft.to_zeta(d)))
                
            d,p,energy = pd.solve_area_pd(warpdata, d, p, iterations, params, params['Lambda'][level], L_normal, imageWidget)
            
    return d



if __name__ == '__main__':
    
    if QtGui.QApplication.instance()==None:
        app=QtGui.QApplication(sys.argv)
        
    plt.close('all')
    
    
    if 'imageWidget' in locals() and imageWidget is not None:
        imageWidget.close()
        imageWidget = None    
        
    npz = np.load('data.npz')
    images = [npz['I1'], npz['I2']]; G = [npz['G1'], npz['G2']]#这个G是一个4x4的矩阵，到底是存储着什么值？
    data = dict(images=images, G=G, K=npz['K'])#内参矩阵K
    d_init = 1.8
    
    I1 = data['images'][0]; I2 = data['images'][1]; K = data['K']
    '''
    G1 = data['G'][0]
    G2 = data['G'][1]
    print(G1,G2)
    sio.savemat('G2.mat',{'G2':G2})
    print('done!')
    time.sleep(9999)
    '''
    d0 = np.ones(I1.shape)*d_init#depth图初始化的值

    minalpha = 0.015             # 3d points must be seen at least under this angle
    maxz = ft.camera_baseline(data['G'][0], data['G'][1])/np.arctan(minalpha/2)    # compute depth value
    #这里计算了depth的最大的深度，但是为什么是计算出来的最大的深度？ 
    params = dict(warps=15,
                  iterations=20,
                  scalefactor=0.75,
                  levels=100,
                  minSize=48,
                  stop_level=0, 
                  check=10,
                  Lambda=0.005,
                  ref=0,
                  epsilon = 0.001,     # huber epsilon
                  minz=0.8,                 # lower threshold for scene depth
                  maxz=maxz)                # maxz: clamp depth to be smaller than this value

                  
    imageWidget = ft.MplWidget()
    imageWidget.show()
    
    plt.figure('I1'); plt.imshow(I1, cmap='gray'); plt.title('Image 1')
    plt.figure('I2'); plt.imshow(I2, cmap='gray'); plt.title('Image 2')
       
    d = run_ctf(params, d0, data, imageWidget)

