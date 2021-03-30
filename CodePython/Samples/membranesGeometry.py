#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 17:03:37 2020

@author: quenot
"""

import numpy as np

from matplotlib import pyplot as plt


def membraneRegularCones(studyDimentions, studyPixelSize):
    
    
    coneSize=int(np.floor(15/studyPixelSize))
    coneRadius=np.sqrt(coneSize**2/2)
    coneMiddle=int(coneSize/2)
    coneHeight=coneSize*1*studyPixelSize*1e-6
    
    Nx=int(coneSize*(np.floor(studyDimentions[0]/coneSize+2)))
    Ny=int(coneSize*(np.floor(studyDimentions[1]/coneSize+2)))
    
    cone=np.zeros((coneSize,coneSize))
    
    for i in range(coneSize):
        for j in range(coneSize):
            #cone[i,j]=1-np.sqrt(abs(0.5-i/coneRadius)**2+abs(0.5-j/coneRadius)**2)
            
            #Speres
#            cone[i,j]=2*np.sqrt(coneRadius**2-(np.sqrt((i-coneMiddle)**2+(j-coneMiddle)**2))**2)
            
            ##Cones
#            cone[i,j]=(coneRadius-(np.sqrt((i-coneMiddle)**2+(j-coneMiddle)**2)))/coneRadius*coneHeight
#            
            
            #Pyramides
            if coneSize-j<=i<j:
                cone[i,j]=(coneSize-j)
            if j<=coneSize-i<=coneSize-j:
                cone[i,j]=j                   
            if i<=j<coneSize-i:
                cone[i,j]=i
            if coneSize-i<=coneSize-j<=i:
                cone[i,j]=(coneSize-i)
                            
            
    cone=cone/coneMiddle*coneHeight
    membrane=np.zeros((Nx,Ny))
 
    for i in range(int(np.floor(Nx/coneSize))):
        for j in range(int(np.floor(Ny/coneSize))):
        
#            membrane[i*coneSize:i*coneSize+coneSize,j*coneSize:j*coneSize+coneSize]=cone
            
            if (i+j)%2==0:
                membrane[i*coneSize:i*coneSize+coneSize,j*coneSize:j*coneSize+coneSize]=coneHeight+cone
   
            if (i+j)%2==1:
                membrane[i*coneSize:i*coneSize+coneSize,j*coneSize:j*coneSize+coneSize]=coneHeight-cone

            

    return membrane[0:studyDimentions[0],0:studyDimentions[1]]+10e-6


if __name__ == "__main__":
    N=600
        
    Membrane=membraneRegularCones([N, N], 2)    
        
    print('\nMembrane')
    plt.figure()
    plt.imshow(imagemedian)
    plt.show()
    
    print("\nProfile")
    plt.figure()
    plt.plot(range(N),Membrane[:,50])
    plt.show()    
    print("\nProfile")
    plt.figure()
    plt.plot(range(N),Membrane[25,:])
    plt.show()