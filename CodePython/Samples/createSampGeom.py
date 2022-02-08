#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 17:37:42 2020

@author: quenot
"""
import numpy as np
from xml.dom import minidom
from InputOutput.pagailleIO import openImage, saveEdf, openSeq
from matplotlib import pyplot as plt
import imutils

def CreateSampleSphere(myName, dimX, dimY, pixelSize):
    """
    Creates simple geometry of a sphere centered on the image

    Args:
        myName (str): name of the sample in the .xml file.
        dimX (int): study dimension in x.
        dimY (int): study dimension in y.
        pixelSize (float): pixel size in the sample plane.

    Returns:
        3D numpy array: thickness maps for different materials (only one in the sphere case) [material, dimX, dimY].
        parameters (dictionnary of tuples): Any parameter related to the geomtry that you want stored in the final text file. (value, "unit")

    """
    parameters={}
    
    #Get additionnal parameters from xml file
    xmlSampleFileName="xmlFiles/Samples.xml"
    xmldocSample = minidom.parse(xmlSampleFileName)
    for currentSample in xmldocSample.documentElement.getElementsByTagName("sample"):
        correctSample = getText(currentSample.getElementsByTagName("name")[0])
        if correctSample == myName:
            myRadius0=float(getText(currentSample.getElementsByTagName("myRadius")[0]))

    if myRadius0/pixelSize*2>max(dimX,dimY):
        print("/!\ Sphere size bigger than the field of view!")
        
    Sample=np.zeros((1,dimX,dimY))
    myRadius=myRadius0/pixelSize
    for i in range(dimX):
        for j in range(dimY):
            dist=(dimX/2-i)**2+(dimY/2-j)**2
            if dist<myRadius**2:
                Sample[0,i,j]=2*np.sqrt(myRadius**2-dist)
                
    # save parameters in the dictionnary
    parameters['Sphere_radius']=(myRadius0, 'um')
    return Sample*pixelSize*1e-6, parameters


def CreateSampleCylindre(myName, dimX, dimY, pixelSize):
    """
    Creates geometry for a cylinder.

    Args:
        myName (str): name of the sample in the .xml file.
        dimX (int): study dimension in x.
        dimY (int): study dimension in y.
        pixelSize (float): pixel size in the sample plane.

    Returns:
        3D numpy array: thickness maps for different materials (only one in the cylinder case) [material, dimX, dimY].
        parameters (dictionnary of tuples): Any parameter related to the geomtry that you want stored in the final text file. (value, "unit")

    """
    parameters={}
    
    #Get additionnal parameters from xml file
    xmlSampleFileName="xmlFiles/Samples.xml"
    xmldocSample = minidom.parse(xmlSampleFileName)
    for currentSample in xmldocSample.documentElement.getElementsByTagName("sample"):
        correctSample = getText(currentSample.getElementsByTagName("name")[0])
        
        if correctSample == myName:
            myRadius_um=float(getText(currentSample.getElementsByTagName("myRadius")[0]))
            Orientation=float(getText(currentSample.getElementsByTagName("myOrientation")[0]))
    
    if myRadius_um/pixelSize*2>max(dimX,dimY):
        print("/!\ Cylinder size bigger than the field of view!")
    
    Sampleb=np.zeros(((1,dimX,dimY)))
    Nxp=2*dimX
    Nyp=2*dimY
    diffx=int((Nxp-dimX)/2)
    diffy=int((Nyp-dimY)/2)

    Sample=np.zeros((Nxp,Nyp))
    myRadius=myRadius_um/pixelSize
    
    for j in range(Nyp):
        if (abs(Nyp/2-j)<myRadius):
            Sample[:,j]=2*np.sqrt(myRadius**2-(Nyp/2-j)**2)
    Samplec=imutils.rotate(Sample, angle=Orientation)
    Sampleb[0]+=Samplec[diffx:diffx+dimX,diffy:diffy+dimY]
    
    #Save parameters in the dictionnary
    parameters['Cylinder_radius']=(myRadius_um, 'um')
    parameters['Cylinder_orientation']=(Orientation, 'degree')
            
    return Sampleb*pixelSize*1e-6, parameters
            

def CreateSampleSpheresInCylinder(myName, dimX, dimY, pixelSize):
    """
    Creates the geometry of a vertical cylinder with two spheres in it.

    Args:
        myName (str): name of the sample in the .xml file.
        dimX (int): study dimension in x.
        dimY (int): study dimension in y.
        pixelSize (float): pixel size in the sample plane.

    Returns:
        3D numpy array: thickness maps for different materials (3 in this case) [material, dimX, dimY].
        parameters (dictionnary of tuples): Any parameter related to the geomtry that you want stored in the final text file. (value, "unit")

    """
    #Define geometry parameters
    myRadius0=500 #um Spheres radius
    myRadius2=myRadius0*2 #um  Cylinder radius
    posY=dimY//2 # horizontal position of the cylinder 
    posXmuscle=int(np.round(myRadius0*3/pixelSize)) #vertical position of 1st sphere
    posXcart=int(np.round(myRadius0*7/pixelSize)) #vertical position of 2nd sphere
    
    parameters={}
    Sample=np.zeros((3,dimX,dimY))
    myRadius=myRadius0/pixelSize
    patchSize2=int(np.ceil(myRadius))
    patchSize=patchSize2*2
    spherePatch=np.zeros((patchSize,patchSize))
    
    for i in range(patchSize):
        for j in range(patchSize):
            dist=(patchSize/2-i)**2+(patchSize/2-j)**2
            if dist<myRadius**2:
                spherePatch[i,j]=2*np.sqrt(myRadius**2-(patchSize/2-j)**2-(patchSize/2-i)**2)
            
    Tube=np.zeros((dimX,dimY))
    myRadius=myRadius2/pixelSize
    for j in range(dimY):
        if (abs(dimY/2-j)<myRadius):
            Tube[:,j]=2*np.sqrt(myRadius**2-(dimY/2-j)**2)
            
    Sample[0,posXmuscle-patchSize2:posXmuscle+patchSize2,posY-patchSize2:posY+patchSize2]=spherePatch
    Sample[1,posXcart-patchSize2:posXcart+patchSize2,posY-patchSize2:posY+patchSize2]=spherePatch
    Sample[2]=Tube-Sample[0]-Sample[1]
    
    plt.figure()
    plt.imshow(Sample[2])
    plt.title('Sample geometry of the cylinder')
    plt.colorbar()
    plt.show()
    
    #Save parameters in the dictionnary
    parameters['Spheres_radius']=(myRadius0, 'um')
    parameters['Cylinder_radius']=(myRadius2, 'um')
    parameters['Position_Sphere_1']=(posXmuscle*pixelSize, 'um')
    parameters['Position_Sphere_2']=(posXcart*pixelSize, 'um')
    return Sample*pixelSize*1e-6, parameters
            
def CreateSampleSpheresInParallelepiped(myName, dimX0, dimY0, pixelSize):  
    """
    Creates the geometry of a vertical parallelepiped with two spheres in it.

    Args:
        myName (str): name of the sample in the .xml file.
        dimX (int): study dimension in x.
        dimY (int): study dimension in y.
        pixelSize (float): pixel size in the sample plane.

    Returns:
        3D numpy array: thickness maps for different materials (3 in this case) [material, dimX, dimY].
        parameters (dictionnary of tuples): Any parameter related to the geomtry that you want stored in the final text file. (value, "unit")

    """
    
    # SAMPLE PARAMETERS
    myRadius0=500 #um radius of inner spheres
    orientation=15 #to tilt the sample
    myRadius2=myRadius0*2 #um width of parallelepipede
    
    
    parameters={}
    margin=max(dimX0, dimY0)//2
    dimX=dimX0+2*margin
    dimY=dimY0+2*margin
    Sample=np.zeros((3,dimX,dimY))
    posY=dimY//2 
    posXmuscle=dimX*2//5 # Positions of the spheres
    posXcart=dimX*3//5
    
    myRadius=myRadius0/pixelSize
    patchSize2=int(np.ceil(myRadius))
    patchSize=patchSize2*2
    spherePatch=np.zeros((patchSize,patchSize))
    
    if abs(posXcart-posXmuscle)<myRadius*2:
        print("/!\ sample spheres overlapping!")
        
    for i in range(patchSize):
        for j in range(patchSize):
            dist=(patchSize/2-i)**2+(patchSize/2-j)**2
            if dist<myRadius**2:
                spherePatch[i,j]=2*np.sqrt(myRadius**2-(patchSize/2-j)**2-(patchSize/2-i)**2)
            
    Tube=np.zeros((dimX,dimY))
    myRadius=myRadius2/pixelSize
    for j in range(dimX):
        if (abs(dimY/2-j)<myRadius*3/4):
            Tube[:,j]=myRadius*2#2*np.sqrt(myRadius**2-(dimX/2-j)**2)
        if myRadius>(j-dimY/2)>=myRadius*3/4:
            Tube[:,j]=myRadius/2*3+2*np.sqrt((myRadius/4)**2-(j-(dimY/2+myRadius*3/4))**2)
        if -myRadius<j-dimY/2<=-myRadius*3/4:
            Tube[:,j]=myRadius/2*3+2*np.sqrt((myRadius/4)**2-(j-(dimY/2-myRadius*3/4))**2)
    
    Sample[0,posXmuscle-patchSize2:posXmuscle+patchSize2,posY-patchSize2:posY+patchSize2]=spherePatch
    Sample[1,posXcart-patchSize2:posXcart+patchSize2,posY-patchSize2:posY+patchSize2]=spherePatch
    Sample[2]=Tube-Sample[0]-Sample[1]
    
    Sample[0]=imutils.rotate(Sample[0], angle=orientation)
    Sample[1]=imutils.rotate(Sample[1], angle=orientation)
    Sample[2]=imutils.rotate(Sample[2], angle=orientation)
    Sample=Sample[:,margin:dimX-margin, margin:dimY-margin]      
    
    plt.figure()
    plt.imshow(Sample[2])
    plt.title("Parallelepipede")
    plt.colorbar()
    plt.show()  
    
    #Save parameters in the dictionnary
    parameters['Spheres_radius']=(myRadius0, 'um')
    parameters['Parallelepipede_size']=(myRadius2, 'um')
    parameters['Position_Sphere_1']=(posXmuscle*pixelSize, 'um')
    parameters['Position_Sphere_2']=(posXcart*pixelSize, 'um')
    return Sample*pixelSize*1e-6, parameters



def CreateYourSampleGeometry(myName, dimX0, dimY0, pixelSize):  
    """
    Creates the geometry EXAMPLE TO MODIFY
    Create one thickness map per material - give as many material in the xml as the thickness maps here

    Args:
        myName (str): name of the sample in the .xml file.
        dimX (int): study dimension in x.
        dimY (int): study dimension in y.
        pixelSize (float): pixel size in the sample plane.

    Returns:
        3D numpy array: thickness maps in m!! for different materials [material, dimX, dimY].
        parameters (dictionnary of tuples): Any parameter related to the geomtry that you want stored in the final text file. (value, "unit")

    """
    thickness=5*1e-6 #5um
    Geometry=np.ones((1,dimX0, dimY0))*thickness
    
    parameters={}
    parameters['geometry thickness']=(thickness, 'um')
    
    return Geometry, parameters

def getText(node):
    return node.childNodes[0].nodeValue