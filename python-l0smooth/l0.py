# -*- coding: utf-8 -*-

import numpy as np
import time
import cv2
import scipy.io as sio
def circshift(new_psf, shift_row, shift_col):
   row, col = new_psf.shape
   shift_row = (row+(shift_row%row))%row
   shift_col = (col+(shift_col%col))%col
   temp = new_psf.copy()   
   if (not shift_row == 0):
      new_psf[0:shift_row, :]=temp[row - shift_row:row, :]
      new_psf[shift_row:row , :]=temp[0:row-shift_row, :]
   if (not shift_col == 0):
      new_psf[:, 0:shift_col ]=temp[:, col - shift_col:col ]
      new_psf[:, shift_col:col ]=temp[:, 0:col-shift_col ]
   return new_psf


def dst2comp(src):
   src_1 = np.zeros((src.shape[0],src.shape[1]), dtype = complex)
   src_1.real = src[:,:,0]
   src_1.imag = src[:,:,1]
   src = np.copy(src_1)
   return src

def psf2otf(psf, outSize):
   '''
   Data:2016-2-16
   Implenment Matlab pst2osf function
   '''
   height,width = psf.shape
   shift_row = -1*(np.floor(height*0.5))
   shift_col = -1*(np.floor(width*0.5))
   new_psf = np.zeros((outSize[0],outSize[1]))
   for i in range(0,height):
      for j in range(0,width):
         new_psf[i][j] = psf[i][j]
         j += 1
      i += 1
   new_psf = circshift(new_psf,shift_row,shift_col) 
   print new_psf;time.sleep(99999)
   otf = np.zeros((outSize[0],outSize[1]), dtype = complex)
   otf1 = cv2.dft(new_psf, flags = cv2.DFT_COMPLEX_OUTPUT)
   otf.real = otf1[:,:,0]
   otf.imag = otf1[:,:,1]
   return otf
   
def l0_norm(Im,lamb=2e-2,kappa=2):
   S = Im*(1/255.0)
   betamax = 1e5
   fx = np.array([[1,-1]])
   fy = np.array([[1],[-1]])
   N,M,D = S.shape
   sizeI2D = np.array([N,M])
   otfFx = psf2otf(fx, sizeI2D)
   
   otfFy = psf2otf(fy, sizeI2D)
   Normin1 = np.zeros((N, M, D),dtype = complex)
   single_channel = np.zeros((D))
   single_channel = cv2.split(S)
   for k in range(3):
      temp_S = single_channel[k] 
      Normin1[:,:,k] = np.fft.fft2(temp_S)
      #Normin1[:,:,k] = cv2.dft(single_channel[k],flags = cv2.DFT_COMPLEX_OUTPUT)
   Denormin2 = abs(otfFx)**2 + abs(otfFy)**2
   beta = 2*lamb
   while beta < betamax:
      Denormin = 1 + beta * Denormin2
#h-v subproblem
      dx = np.zeros((N, M, D))
      dy = np.zeros((N, M, D))
      shifted_x = np.zeros((N, M))      
      shifted_y = np.zeros((N, M))      
   
      for k in range(D):
         shifted_x = np.copy(single_channel[k]) 
         shifted_x = circshift(shifted_x, 0, -1)
         dx[:,:,k] = shifted_x - single_channel[k]
         shifted_y = np.copy(single_channel[k])
         shifted_y = circshift(shifted_y, -1, 0)
         dy[:,:,k] = shifted_y - single_channel[k]
      
      for i in range(N):
         for j in range(M):
            val = (dx[i,j,0])**2+(dy[i,j,0])**2+(dx[i,j,1])**2+(dy[i,j,1])**2+(dx[i,j,2])**2+(dy[i,j,2])**2
            if (val <lamb/beta):
               dx[i,j,0] = dx[i,j,1] = dx[i,j,2] = 0.0
               dy[i,j,0] = dy[i,j,1] = dy[i,j,2] = 0.0
#S subproblem

      for k in range(D):
         shift_dx = np.copy(dx[:,:,k])
         shift_dx = circshift(shift_dx, 0 ,1)
         ddx = shift_dx - dx[:,:,k]
      
         shift_dy = np.copy(dy[:,:,k])
         shift_dy = circshift(shift_dy, 1, 0)
         ddy = shift_dy - dy[:,:,k]
         
         Normin2 = ddx + ddy
         FNormin2 = cv2.dft(Normin2, flags = cv2.DFT_COMPLEX_OUTPUT)
         FNormin2 = dst2comp(FNormin2)
         FS = Normin1[:,:,k] + beta * FNormin2
         #FS = np.zeros((N, M, D),dtype = complex)
         for i in range(N):
            for j in range(M):
               FS[i,j] = FS[i,j]/Denormin[i,j]
         FS_comp = np.zeros((N,M,2))
         FS_comp[:,:,0] = FS.real
         FS_comp[:,:,1] = FS.imag
         ifft = cv2.idft(FS_comp, flags = cv2.DFT_SCALE |cv2.DFT_COMPLEX_OUTPUT)
         #ifft = np.fft.ifft(FS)
         single_channel[k] = ifft[:,:,0]

      beta = beta * kappa
      #beta = betamax
      print('.')
      S = cv2.merge(single_channel,3)
   
   return S

if __name__ == '__main__':
   '''
   img= cv2.imread('pflower.jpg')
   cv2.namedWindow('Image')
   cv2.imshow('Image',img)
   cv2.waitKey(0)
   cv2.destroyAllWindows() 
   print(img.shape)
   '''   
   Im = cv2.imread('pflower.jpg')
   image = l0_norm(Im)
   cv2.namedWindow('L0_norm')
   cv2.imshow('L0_norm',image)
   cv2.waitKey(0)
   cv2.destoryAllWindows()
    
