#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 16:47:58 2021

@author: quenot
"""
import numpy as np
from numba import jit 
from scipy.ndimage import gaussian_filter
from matplotlib import pyplot as plt


def gaussian_shape(sigma):
    dim=round(sigma*3)*2+1
    
    Qx, Qy = np.meshgrid((np.arange(0, dim) - np.floor(dim / 2) ), (np.arange(0, dim) - np.floor(dim / 2)) ) #frequency ranges of the images in fqcy space

    #sigmaY = sig_scale

    g = np.exp(-(((Qx)**2) / 2. / sigma**2 + ((Qy)**2) / 2. / sigma**2))

    return g/np.sum(g)

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
    margin2=15
    # print("Sum insidentInsentity beginning function refraction", np.sum(intensityRefracted))
    
    #Get from the phase information to the displacement in pixel after propagation
    dphix,dphiy=np.gradient(phi,studyPixelSize*1e-6, edge_order=2)
    Dx=dphix*propagationDistance/k/(studyPixelSize*1e-6*magnification)
    Dy=dphiy*propagationDistance/k/(studyPixelSize*1e-6*magnification)
    
    #Get rid of aberrant values or values that wont influence the result
    Dx[abs(Dx)<1e-12]=0
    Dy[abs(Dy)<1e-12]=0
    intensityRefracted[abs(Dx)>Nx]=0
    intensityRefracted[abs(Dy)>Ny]=0
    Dx[abs(Dx)>Nx]=0
    Dy[abs(Dy)>Ny]=0
    Dx=np.pad(Dx,margin2,mode='constant')
    Dy=np.pad(Dy,margin2,mode='constant')
    intensityRefracted=np.pad(intensityRefracted,margin2,mode='constant')
    
    #initialize the resulting intensity matrix with a margin for the calculation
    intensityRefracted2=np.zeros((Nx+2*margin2, Ny+2*margin2))
    
    DxFloor=Dx.astype(np.int)
    DyFloor=Dy.astype(np.int)
    
    
    #Call the fast function for the loop on all the pixels
    intensityRefracted2=fastloopNumba(Nx+2*margin2, Ny+2*margin2,intensityRefracted,intensityRefracted2,Dy,Dx,DxFloor, DyFloor)
    intensityRefracted2 = intensityRefracted2[margin2:Nx+margin2,margin2:Ny+margin2]
    
    #check for errors in the calculation (there shouldn't be but we never know)
    if np.isnan(intensityRefracted2).any() or np.any((abs(intensityRefracted2) > 1e50)):
        raise Exception("The calculated intensity refractive includes some nans or insane values")
    
    # print("Sum insidentInsentity end function refraction", np.sum(intensityRefracted2))
    
    return intensityRefracted2,Dx, Dy   

def fastRefractionDF(intensityRefracted, phi, propagationDistance,Energy, magnification,studyPixelSize, darkField):
    """
    Calculates the intensity after propagation from the angles of refraction

    Args:
        intensityRefracted (2D numpy array): Intensity before propagation.
        phi (2D numpy array): phase information.
        propagationDistance (Float): propagation distance (in m).
        Energy (Float): energy of the spectrum currently computed (in keV).
        magnification (Float): magnification from the source to the distance before porpagation and after.
        studyPixelSize (Float): voxel size on the plane before propagation (in um).
        darkField (2D numpy array): averaage scattering angle (rad)

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
    darkField=darkField*propagationDistance/(studyPixelSize*1e-6*magnification)
    maxDF=np.max(darkField)
    # print('Max df in pixels:', maxDF)
    margin2=int(np.ceil(maxDF*6))
    
    #Get from the phase information to the displacement in pixel after propagation
    dphix,dphiy=np.gradient(phi,studyPixelSize*1e-6, edge_order=2)
    Dx=dphix*propagationDistance/k/(studyPixelSize*1e-6*magnification)
    Dy=dphiy*propagationDistance/k/(studyPixelSize*1e-6*magnification)
    
    #Get rid of aberrant values or values that wont influence the result
    Dx[abs(Dx)<1e-12]=0
    Dy[abs(Dy)<1e-12]=0
    # darkField[abs(darkField)<1e-12]=0
    intensityRefracted[abs(Dx)>Nx]=0
    intensityRefracted[abs(Dy)>Ny]=0
    Dx[abs(Dx)>Nx]=0
    Dy[abs(Dy)>Ny]=0
    Dx=np.pad(Dx,margin2,mode='constant')
    Dy=np.pad(Dy,margin2,mode='constant')
    darkField[darkField>Nx/4]=0
    
    darkField=np.pad(darkField,margin2,mode='constant')
    intensityRefracted=np.pad(intensityRefracted,margin2,mode='constant')

    #initialize the resulting intensity matrix with a margin for the calculation
    intensityRefracted2=np.zeros((Nx+2*margin2, Ny+2*margin2))
    intensityRefracted2DF=np.zeros((Nx+2*margin2, Ny+2*margin2))

    DxFloor=Dx.astype(np.int)
    DyFloor=Dy.astype(np.int)
    
    #Call the fast function for the loop on all the pixels
    intensityNoDF=np.copy(intensityRefracted)
    intensityNoDF[darkField!=0]=0
    intensityDF=np.copy(intensityRefracted)
    intensityDF[darkField==0]=0

    plt.figure()
    plt.imshow(darkField)
    plt.title('dark-field')
    plt.colorbar()
    plt.show()
    
    intensityRefracted2=fastloopNumba(Nx+2*margin2, Ny+2*margin2,intensityNoDF,intensityRefracted2,Dy,Dx,DxFloor, DyFloor)
    intensityRefracted2DF=fastloopNumba(Nx+2*margin2, Ny+2*margin2,intensityDF,intensityRefracted2DF,Dy,Dx,DxFloor, DyFloor)
    # darkField=fastloopNumba(Nx+2*margin2, Nx+2*margin2,darkField,darkField,Dy,Dx,DxFloor, DyFloor)
    intensityRefracted3=np.zeros((Nx+2*margin2, Ny+2*margin2))
    
    plt.figure()
    plt.imshow(darkField)
    plt.title('dark-field')
    plt.colorbar()
    plt.show()
    
    # print('max intensityRefracted', np.max(intensityRefracted))
    
    for i in range(margin2,Nx+margin2):
        for j in range(margin2,Ny+margin2):
            if intensityRefracted2DF[i,j]!=0:
                if darkField[i,j]!=0:
                    # print(i,j)
                    currDF=darkField[i,j]/2#/0.67
                    patch=gaussian_shape(currDF)
                    size2=patch.shape[0]//2
                    patch=patch*intensityRefracted2DF[i,j]
                    # print("Max patch, DF, Dx, Dy", np.max(patch), currentDF,Dx[i,j],Dy[i,j] )
                    #print(i,j, currDF, size2,patch.shape[0])
                    intensityRefracted3[i-size2:i+size2+1,j-size2:j+size2+1]+=patch
                else:
                    intensityRefracted3[i,j]+=intensityRefracted2DF[i,j]
    
    intensityRefracted3+=intensityRefracted2
    intensityRefracted3 = intensityRefracted3[margin2:Nx+margin2,margin2:Ny+margin2]
    
    #check for errors in the calculation (there shouldn't be but we never know)
    if np.any((abs(intensityRefracted3) > 1e50)):
        raise Exception("The calculated intensity refractive includes some  insane values")
    if np.isnan(intensityRefracted3).any() :
        raise Exception("The calculated intensity refractive includes some nans")
    
    # print("Sum insidentInsentity end function refraction", np.sum(intensityRefracted3))
    return intensityRefracted3,Dx, Dy  

@jit(nopython=True)
def fastloopNumba(Nx, Ny,intensityRefracted,intensityRefracted2,Dy,Dx,DxFloor, DyFloor):
    """
    Accelerated part of the refraction calculation

    Args:
        Nx (int): shape[0] of the intensity refracted.
        Ny (int): shape[1] of the intensity refracted.
        intensityRefracted (2d numpy array): intensity before propag.
        intensityRefracted2 (2d numpy array): intensity after propag.
        Dy (2d numpy array): Displacement along x (in voxel).
        Dx (2d numpy array): Displacement along y (in voxel).
        DxFloor (2d numpy array): floored displacement.
        DyFloor (2d numpy array): floored displacement.

    Returns:
        intensityRefracted2 (2d numpy array): intensity after propag.

    """
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


@jit(nopython=True)
def fastloopNumbaDF(Nx, Ny,intensityRefracted,intensityRefracted2,Dy,Dx,DxFloor, DyFloor, DF):
    """
    Accelerated part of the refraction calculation

    Args:
        Nx (int): shape[0] of the intensity refracted.
        Ny (int): shape[1] of the intensity refracted.
        intensityRefracted (2d numpy array): intensity before propag.
        intensityRefracted2 (2d numpy array): intensity after propag.
        Dy (2d numpy array): Displacement along x (in voxel).
        Dx (2d numpy array): Displacement along y (in voxel).
        DxFloor (2d numpy array): floored displacement.
        DyFloor (2d numpy array): floored displacement.

    Returns:
        intensityRefracted2 (2d numpy array): intensity after propag.

    """
    for i in range(Nx):
        for j in range(Ny):
            Iij=intensityRefracted[i,j]
            Dxtmp=Dx[i,j]
            Dytmp=Dy[i,j]
            if DF[i,j]!=0:
                continue
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
