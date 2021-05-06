#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 16:47:58 2021

@author: quenot
"""
import numpy as np
from ctypes import *
# fastloopm=CDLL("fastloopmodule.so")
# from fastloop import fastloop
from time import time
import ctypes
from scipy import ndimage, LowLevelCallable
from numba import jit 

def fastRefraction(intensityRefracted, phi, propagationDistance,Energy, magnification,studyPixelSize):
    """
    Calculates the intensity after propagation from the angles of refraction

    Args:
        intensityRefracted (2D numpy array): Intensity before propagation.
        phi (2D numpy array): phase information.
        propagationDistance (Float): propagation distance (in m).
        Energy (Float): energy of the spectrum currently computed (in keV).
        magnification (Float): magnification from the source to the distance before porpagation and after.
        studyPixelSize (Float): voxel size on the plane before propagation (in um).

    Raises:
        Exception: If calculations give aberrant values.

    Returns:
        intensityRefracted2 (2D numpy array): intensity after propagation.
        Dx (2D numpy array): Real displacement map calculated.
        Dy (2D numpy array): Real displacement map calculated.

    """
    #Get usefull experimental constants
    Lambda=6.626*1e-34*2.998e8/(Energy*1000*1.6e-19)
    k=2*np.pi/Lambda
    Nx,Ny=intensityRefracted.shape
    margin2=2
    
    #Get from the phase information to the displacement in pixel after propagation
    dphix,dphiy=np.gradient(phi,studyPixelSize*1e-6, edge_order=2)
    Dx=dphix*propagationDistance/k/(studyPixelSize*1e-6*magnification)
    Dy=dphiy*propagationDistance/k/(studyPixelSize*1e-6*magnification)
    
    #Get rid of aberrant values or values that wont influence the result
    Dx[abs(Dx)<1e-12]=0
    Dy[abs(Dy)<1e-12]=0
    intensityRefracted[abs(Dx)>1e3]=0
    intensityRefracted[abs(Dy)>1e3]=0
    Dx[abs(Dx)>1e3]=0
    Dy[abs(Dy)>1e3]=0
    
    #initialize the resulting intensity matrix with a margin for the calculation
    intensityRefracted2=np.zeros((Nx+2*margin2, Ny+2*margin2))
    
    DxFloor=Dx.astype(np.int)
    DyFloor=Dy.astype(np.int)
    
    #Call the fast function for the loop on all the pixels
    fastloopNumba(Nx, Ny,intensityRefracted,intensityRefracted2,Dy,Dx,DxFloor, DyFloor)
    intensityRefracted2 = intensityRefracted2[margin2:Nx+margin2,margin2:Ny+margin2]
    
    #check for errors in the calculation (there shouldn't be but we never know)
    if np.isnan(intensityRefracted2).any() or np.any((abs(intensityRefracted2) > 1e50)):
        raise Exception("The calculated intensity refractive includes some nans or insane values")
    
    return intensityRefracted2,Dx, Dy    

@jit(nopython=True)
def fastloopNumba(Nx, Ny,intensityRefracted,intensityRefracted2,Dy,Dx,DxFloor, DyFloor):
    for i in range(Nx):
        for j in range(Ny):
            Iij=intensityRefracted[i,j]
            Dxtmp=Dx[i,j]
            Dytmp=Dy[i,j]
            if Dxtmp==0 and Dytmp==0:
                intensityRefracted2[i,j]+=Iij
                continue
            inew=i
            jnew=j
            #Calculating displacement bigger than a pixel
            if abs(Dxtmp)>1:
                inew=i+int(np.floor(Dxtmp))
                Dxtmp=Dxtmp-np.floor(Dxtmp)
            if abs(Dytmp)>1:
                jnew=j+int(np.floor(Dytmp))
                Dytmp=Dytmp-np.floor(Dytmp)
            #Calculating sub-pixel displacement
            if 0<=inew<Nx:
                if 0<=jnew<Ny:
                    intensityRefracted2[inew,jnew]+=Iij*(1-abs(Dxtmp))*(1-abs(Dytmp))
                    if inew<Nx-1:
                        if Dxtmp>=0:
                            if jnew<Ny-1:
                                if Dytmp>=0:
                                    intensityRefracted2[inew+1,jnew]+=Iij*Dxtmp*(1-Dytmp)
                                    intensityRefracted2[inew+1,jnew+1]+=Iij*Dxtmp*Dytmp
                                    intensityRefracted2[inew,jnew+1]+=Iij*(1-Dxtmp)*Dytmp
                            if jnew>0:
                                if Dytmp<0:
                                    intensityRefracted2[inew+1,jnew]+=Iij*Dxtmp*(1-abs(Dytmp))
                                    intensityRefracted2[inew+1,jnew-1]+=Iij*Dxtmp*abs(Dytmp)
                                    intensityRefracted2[inew,jnew-1]+=Iij*(1-Dxtmp)*abs(Dytmp)
                    
                    if inew>0:
                        if Dxtmp<0:
                            if jnew<Ny-1:
                                if Dytmp>=0:
                                    intensityRefracted2[inew-1,jnew]+=Iij*abs(Dxtmp)*(1-Dytmp)
                                    intensityRefracted2[inew-1,jnew+1]+=Iij*abs(Dxtmp)*Dytmp
                                    intensityRefracted2[inew,jnew+1]+=Iij*(1-abs(Dxtmp))*Dytmp
                            if jnew>0:
                                if Dytmp<0:
                                    intensityRefracted2[inew-1,jnew]+=Iij*abs(Dxtmp)*(1-abs(Dytmp))
                                    intensityRefracted2[inew-1,jnew-1]+=Iij*Dxtmp*Dytmp
                                    intensityRefracted2[inew,jnew-1]+=Iij*(1-abs(Dxtmp))*abs(Dytmp)
    return intensityRefracted2


