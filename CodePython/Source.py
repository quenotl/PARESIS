#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 17:30:05 2020

@author: quenot
"""

from xml.dom import minidom
import numpy as np
import matplotlib.pyplot as plt
from spekpy import Spek


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
        self.myVoltage=0

    def __str__(self):
        return f'Source: {self.myName}\n Size: {self.mySize} um\n Type: {self.myType}\n Spectrum: {self.mySpectrum} (KeV, 1) \n'
        
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
    
    def totalFlux(self):
        return sum([flux for _, flux in self.mySpectrum])
            
    def setMySpectrum(self):
        """
        sets the source spectrum from xml file value for monochromatic or Spekpy for polychromatic

        Returns:
            None.

        """
        # Monochromatic source case
        if self.myType=="Monochromatic":
            self.mySpectrum.append((float(self.getText(self.currentSource.getElementsByTagName("myEnergy")[0])),1))
            return
        
        # Polychromatic source case
        if self.myType=="Polychromatic":
            #get spectrum from spekpy
            s = Spek(kvp=self.myVoltage, targ=self.myTargetMaterial, dk=self.myEnergySampling)
            #taking into account exiting filter?
            if self.exitingWindowMaterial is not None:
                s.filter(self.exitingWindowMaterial, self.exitingWindowThickness)
            self.mySpectrum = list(zip(*s.get_spectrum()))
            plt.figure('Source Filtered Spectrum')
            plt.plot(*list(zip(*self.mySpectrum)))
            plt.xlabel('Energy (keV)')
            plt.title("Source filtered spectrum")
            plt.show(block=False)
            
            while self.mySpectrum[0][1]/s.get_flu()<0.001:
                self.mySpectrum.pop(0)

            plt.figure('Resampled Spectrum')
            plt.plot(*list(zip(*self.mySpectrum)))
            plt.xlabel('Energy (keV)')
            plt.title("Resampled spectrum")
            plt.show(block=True)
            return
            
        raise ValueError("type of source not recognized")
        
    
    def getText(self,node):
        return node.childNodes[0].nodeValue
    
    
if __name__ == "__main__":
    s = Spek(kvp=40,th=12, targ='W')
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
    