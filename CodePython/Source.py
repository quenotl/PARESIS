#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 17:30:05 2020

@author: quenot
"""

import os
import glob
from xml.dom import minidom
import numpy as np
from matplotlib import pyplot as plt
from scipy.ndimage.filters import gaussian_filter, median_filter
import spekpy as sp



class Source:
    def __init__(self):
        self.xmlSourcesFileName="xmlFiles/Sources.xml"
        self.xmldocSources = minidom.parse(self.xmlSourcesFileName)
        self.myName=""
        self.mySpectrumPath=""
        self.mySpectrum=[]
        self.mySize=0.
        self.myType=None
        
    def defineCorrectValuesSource(self):
        for currentSource in self.xmldocSources.documentElement.getElementsByTagName("source"):
            correctSource = self.getText(currentSource.getElementsByTagName("name")[0])
            if correctSource == self.myName:
                self.currentSource=currentSource
                self.mySize=float(self.getText(currentSource.getElementsByTagName("mySize")[0]))
                self.myType=self.getText(currentSource.getElementsByTagName("myType")[0])
                if self.myType=="Polychromatic":
                    self.myEnergySampling=float(self.getText(currentSource.getElementsByTagName("myEnergySampling")[0]))
                    self.myVoltage=float(self.getText(currentSource.getElementsByTagName("sourceVoltage")[0]))
                if self.myType=="Monochromatic":
                    self.myEnergySampling=1
                return
            
        raise ValueError("Source not found in the xml file")
            
    def setMySpectrum(self):
        
#        print("type de source:", self.myType)
        if self.myType=="Monochromatic":
            
            self.mySpectrum.append((float(self.getText(self.currentSource.getElementsByTagName("myEnergy")[0])),1))
            spectrum=self.mySpectrum
            return
        
        if self.myType=="Polychromatic":
            
            s = sp.Spek(kvp=self.myVoltage,th=12)
            s.filter('Be', 0.2)
            spectrum=s.get_spectrum()
            
            plt.figure()
            plt.plot(spectrum[0],spectrum[1])
            plt.xlabel('Energy (keV)')
            plt.show()
            
            energyplot=[]
            weightplot=[]
            
            Nen=len(spectrum[0])
            Nbin=int(np.ceil(Nen/self.myEnergySampling/2))
            n=0
            for i in range(Nbin-1):
                currBin=0
                weightBin=0
                energyBin=0
                while currBin<self.myEnergySampling:
                    weightBin+=spectrum[1][n]
                    energyBin+=spectrum[1][n]*spectrum[0][n]
                    n+=1
                    currBin=currBin+0.5
                self.mySpectrum.append((energyBin/weightBin,weightBin))
                energyplot.append(energyBin/weightBin)
                weightplot.append(weightBin)
                
            currBin=0
            weightBin=0
            energyBin=0
            while n<Nen:
                weightBin+=spectrum[1][n]
                energyBin+=spectrum[1][n]*spectrum[0][n]
                n+=1
            self.mySpectrum.append((energyBin/weightBin,weightBin))
            energyplot.append(energyBin/weightBin)
            weightplot.append(weightBin)
                
            plt.figure()
            plt.plot(energyplot,weightplot)
            plt.xlabel('Energy (keV)')
            plt.show()
            
            return
            
        raise ValueError("type of source not recognized")
        
    
    def getText(self,node):
        return node.childNodes[0].nodeValue
    