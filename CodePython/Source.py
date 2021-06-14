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
        self.mySpectrum=[]
        self.mySize=0.
        self.myType=None
        self.exitingWindowMaterial=None
        self.myTargetMaterial='W'
        self.myEnergySampling=1
        
    def defineCorrectValuesSource(self):
        """
        gets all the source parameters from the xml file

        Raises:
            ValueError: Source not found in the xml file.

        Returns:
            None.

        """
        for currentSource in self.xmldocSources.documentElement.getElementsByTagName("source"):
            correctSource = self.getText(currentSource.getElementsByTagName("name")[0])
            if correctSource == self.myName:
                self.currentSource=currentSource
                self.mySize=float(self.getText(currentSource.getElementsByTagName("mySize")[0]))
                self.myType=self.getText(currentSource.getElementsByTagName("myType")[0])
                if self.myType=="Polychromatic":
                    self.myEnergySampling=float(self.getText(currentSource.getElementsByTagName("myEnergySampling")[0]))
                    self.myVoltage=float(self.getText(currentSource.getElementsByTagName("sourceVoltage")[0]))
                    for node in currentSource.childNodes:
                        if node.localName=="exitingWindowMaterial":
                            self.exitingWindowMaterial=self.getText(currentSource.getElementsByTagName("exitingWindowMaterial")[0])
                            self.exitingWindowThickness=float(self.getText(currentSource.getElementsByTagName("exitingWindowThickness")[0]))                            
                        if node.localName=="myTargetMaterial": #Option to chose the target material
                            self.myTargetMaterial=self.getText(currentSource.getElementsByTagName("myTargetMaterial")[0])
                if self.myType=="Monochromatic":
                    self.myEnergySampling=1
                return
            
        raise ValueError("Source not found in the xml file")
            
    def setMySpectrum(self):
        """
        sets the source spectrum from xml file value for monochromatic or Spekpy for polychromatic

        Returns:
            None.

        """
#        print("type de source:", self.myType)

        # Monochromatic source case
        if self.myType=="Monochromatic":
            self.mySpectrum.append((float(self.getText(self.currentSource.getElementsByTagName("myEnergy")[0])),1))
            spectrum=self.mySpectrum
            return
        
        # Polychromatic source case
        if self.myType=="Polychromatic":
            #get spectrum from spekpy
            s = sp.Spek(kvp=self.myVoltage,th=12, targ=self.myTargetMaterial)
            #taking into account exiting filter?
            if self.exitingWindowMaterial is not None:
                s.filter(self.exitingWindowMaterial, self.exitingWindowThickness)
            spectrum=s.get_spectrum()
            
            plt.figure()
            plt.plot(spectrum[0],spectrum[1])
            plt.xlabel('Energy (keV)')
            plt.title("Source filtered spectrum")
            plt.show()
            
            #re-sampling at sourcre energy sampling (to get faster calculation at the end)
            energyplot=[]
            weightplot=[]
            Nen=len(spectrum[0])
            Nbin=int(np.ceil(Nen/self.myEnergySampling/2))
            n=0
            totWeight=0
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
                totWeight+=weightBin
                
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
            
            
            k=0
            while self.mySpectrum[0][1]/totWeight<0.001:
                self.mySpectrum.pop(0)
                energyplot.pop(0)
                weightplot.pop(0)
                k+=1
                
            plt.figure()
            plt.plot(energyplot,weightplot)
            plt.xlabel('Energy (keV)')
            plt.title("Resampled spectrum")
            plt.show()
            
            return
            
        raise ValueError("type of source not recognized")
        
    
    def getText(self,node):
        return node.childNodes[0].nodeValue
    
    
if __name__ == "__main__":
    s = sp.Spek(kvp=40,th=12, targ='W')
    s.filter('Be', 0.2)
    spectrum=s.get_spectrum()
    
    sumEn=0
    sumFlux=0
    stdDev=0
    
    for i in range(len(spectrum[0])):
        sumEn+=spectrum[0][i]*spectrum[1][i]
        sumFlux+=spectrum[1][i]
    meanEn=sumEn/sumFlux
    for i in range(len(spectrum[0])):
        stdDev+=(meanEn-spectrum[0][i])**2*spectrum[1][i]
    
    stdDev=np.sqrt(stdDev/(sumFlux))
    
    print("Energy mean:", sumEn/sumFlux)
    print("Energy stdDev:", stdDev)
    
    plt.figure()
    plt.plot(spectrum[0],spectrum[1]/np.max(spectrum[1]))
    plt.xlabel('Energy (keV)')
    # plt.savefig('/Users/quenot/Library/Mobile Documents/com~apple~CloudDocs/Thèse/ComparaisonSimu/W40kVp.png')
    plt.show()
    
    
    
    
    
