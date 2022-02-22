#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 17:37:42 2020

@author: quenot
"""
import numpy as np
from xml.dom import minidom
import imutils
import matplotlib.pyplot as plt


def CreateSampleCylindre(myName, dimX, dimY, pixelSize):

    xmlSampleFileName = "xmlFiles/Samples.xml"
    xmldocSample = minidom.parse(xmlSampleFileName)

    for currentSample in xmldocSample.documentElement.getElementsByTagName("sample"):
        correctSample = getText(currentSample.getElementsByTagName("name")[0])

        if correctSample == myName:
            myRadius = float(
                getText(currentSample.getElementsByTagName("myRadius")[0]))
            WireAngle = float(
                getText(currentSample.getElementsByTagName("myOrientation")[0]))

    Sampleb = np.zeros(((1, dimX, dimY)))
    Nxp = 2*dimX
    Nyp = 2*dimY
    diffx = int((Nxp-dimX)/2)
    diffy = int((Nyp-dimY)/2)

    Sample = np.zeros((Nxp, Nyp))
    myRadius = myRadius/pixelSize

    for j in range(Nyp):
        if (abs(Nyp/2-j) < myRadius):
            Sample[:, j] = 2*np.sqrt(myRadius**2-(Nyp/2-j)**2)

    Samplec = imutils.rotate(Sample, angle=WireAngle)
    # print("Fylon Wire Geometry")
    plt.figure('Samplec')
    plt.imshow(Samplec)
    plt.colorbar()
    plt.show(block=False)

    Sampleb[0] += Samplec[diffx:diffx+dimX, diffy:diffy+dimY]

    return Sampleb*pixelSize*1e-6, myRadius


def getText(node):
    return node.childNodes[0].nodeValue


if __name__ == '__main__':
    Sample = CreateSampleCylindre('filNylon', 200, 200, 3)
