#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 11:15:32 2021

@author: quenot
"""
import cv2
import glob
from InputOutput.pagailleIO import saveEdf,openImage
import os
import numpy as np

oversamp=4
Ny=int(np.floor(1500*4))
Nx=int(np.floor(500*4))
angle =90
pixelSize=100/oversamp #um
#geometry=generateContrastPhantom(Nx, Ny, pixelSize, angle)

pixsizeStr='%4.1f'%(pixelSize*oversamp)
filepath='/Users/quenot/Documents/Simulations/CodePython_18012021/Samples/ContrastPhantom/CP_21p4_angle90/'
#filepath='/Users/quenot/Documents/Simulations/CodePython_18012021/Samples/Membranes/membraneFe90_21p4um_overSampx10.edf'
imPaths= glob.glob(filepath+'/*') 
imPaths.sort()
Savefilepath="ContrastPhantom/CP_Simap_"+str(pixelSize*oversamp).replace(".","p")+"overSampx"+str(oversamp)+"_angle"+str(angle)
#Savefilepath='/Users/quenot/Documents/Simulations/CodePython_18012021/Samples/Membranes/membraneFe90_21p4um_overSampx'+str(oversamp)
#os.mkdir(Savefilepath)
#geometry=cv2.resize(openImage(filepath),(Ny,Nx))
#saveEdf(geometry,Savefilepath+'.edf')

#if os.path.isfile(filepath+"_mat0"):
for i in range(13):
    filepath=imPaths[i]
    num='%02d'%i
    geometry=cv2.resize(openImage(filepath),(Ny,Nx))
    saveEdf(geometry,Savefilepath+"/mat"+num+".edf")
        
#else:
#     raise Exception("The contrast phantom geometry you are trying to load does not exist or is incorrectly named", filepath)   

oversamp=4
pixelSize=50/oversamp #um
Nx=int(np.floor(5428*0.4))
Ny=int(np.floor(17166*0.4))
pixsizeStr='%4.1f'%(pixelSize*oversamp)
Savefilepath='/Users/quenot/Documents/Simulations/CodePython_18012021/Samples/Membranes/membraneFe90_Simap_21p4um_overSampx'+str(oversamp)
filepath='/Users/quenot/Documents/Simulations/CodePython_18012021/Samples/Membranes/membraneFe90_21p4um_overSampx'+str(oversamp)+'.edf'
#os.mkdir(Savefilepath)
geometry=cv2.resize(openImage(filepath),(Ny,Nx))
saveEdf(geometry,Savefilepath+'.edf')
