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
import xlrd

class Source:
    def __init__(self):
        self.xmlSourcesFileName="xmlFiles/Sources.xml"
        self.xmldocSources = minidom.parse(self.xmlSourcesFileName)
        self.myName=""
        self.mySpectrum=[]
        self.source_dict={}
        self.source_dict["mySize"]=0.
        self.source_dict["myEnergySampling"]=1
        self.source_dict["myType"]=None
        self.spectrumFromXls=False
        
        #initialize units
        self.source_dict["mySize_unit"]="um"
        self.source_dict["myEnergySampling_unit"]="keV"
        self.source_dict["myVoltage_unit"]="kVp"
        self.source_dict["Energy_unit"]="keV"
        self.source_dict["filterThickness_unit"]="mm"
        
        
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
                self.source_dict["mySize"]=float(self.getText(currentSource.getElementsByTagName("mySize")[0]))
                self.source_dict["myType"]=self.getText(currentSource.getElementsByTagName("myType")[0])
                if self.source_dict["myType"]=="Polychromatic":
                    self.source_dict["filterMaterial"]=None
                    self.source_dict["myEnergySampling"]=float(self.getText(currentSource.getElementsByTagName("myEnergySampling")[0]))
                    for node in currentSource.childNodes:
                        if node.localName=="sourceVoltage":
                            self.source_dict["myVoltage"]=float(self.getText(currentSource.getElementsByTagName("sourceVoltage")[0]))
                        if node.localName=="spectrumFromXls":
                            self.spectrumFromXls=bool(self.getText(currentSource.getElementsByTagName("spectrumFromXls")[0]))
                            self.source_dict["pathXlsSpectrum"]=self.getText(currentSource.getElementsByTagName("pathXlsSpectrum")[0])
                            self.source_dict["energyUnit"]=self.getText(currentSource.getElementsByTagName("energyUnit")[0])
                            self.source_dict["energyColumnKey"]=self.getText(currentSource.getElementsByTagName("energyColumnKey")[0])
                            self.source_dict["fluenceColumnKey"]=self.getText(currentSource.getElementsByTagName("fluenceColumnKey")[0])
                        if node.localName=="filterMaterial":
                            self.source_dict["filterMaterial"]=self.getText(currentSource.getElementsByTagName("filterMaterial")[0])
                            self.source_dict["filterThickness"]=float(self.getText(currentSource.getElementsByTagName("filterThickness")[0]))                            
                        if node.localName=="myTargetMaterial": #Option to chose the target material
                            self.source_dict["myTargetMaterial"]=self.getText(currentSource.getElementsByTagName("myTargetMaterial")[0])
                if self.source_dict["myType"]=="Monochromatic":
                    self.source_dict["myEnergySampling"]=1
                    self.source_dict["Energy"]=float(self.getText(self.currentSource.getElementsByTagName("myEnergy")[0]))
                return
            
        raise ValueError("Source not found in the xml file")
            
    def setMySpectrum(self, flu_fluEn=True):
        """
        sets the source spectrum from xml file value for monochromatic or Spekpy for polychromatic

        Returns:
            None.

        """
#        print("type de source:", self.source_dict["myType"])

        # Monochromatic source case
        if self.source_dict["myType"]=="Monochromatic":
            self.mySpectrum.append((self.source_dict["Energy"],1))
            spectrum=self.mySpectrum
            return
        
        # Polychromatic source case
        if self.source_dict["myType"]=="Polychromatic":
            if not self.spectrumFromXls:
                if not "myTargetMaterial" in self.source_dict:
                    self.source_dict["myTargetMaterial"]='W'
                # Polychromatic source case
                #get spectrum from spekpy
                s = sp.Spek(kvp=self.source_dict["myVoltage"],th=12, targ=self.source_dict["myTargetMaterial"],dk=self.source_dict["myEnergySampling"])
                #taking into account exiting filter?
                if self.source_dict["filterMaterial"] is not None:
                    s.filter(self.source_dict["filterMaterial"], self.source_dict["filterThickness"])
                spectrum=s.get_spectrum( flu=flu_fluEn)
                # print(f'mean energy emmited by source: {s.get_emean()}')
                fl=0
                for i in range(len(spectrum[0])):
                    if np.isnan(spectrum[1][i]):
                        spectrum[1][i]=0
                    fl+=spectrum[1][i]
                    
                    
                energyplot=[]
                weightplot=[]
                sumW=0
                for i in range(len(spectrum[0])):
                    if spectrum[1][i]/fl>0.0001:
                        self.mySpectrum.append((spectrum[0][i],spectrum[1][i]/fl))#*self.source_dict["myEnergySampling"]))
                        energyplot.append(spectrum[0][i])
                        weightplot.append(spectrum[1][i]/fl)#*self.source_dict["myEnergySampling"])
                        sumW+=spectrum[1][i]/fl#*self.source_dict["myEnergySampling"]
                
                plt.figure()
                plt.plot(energyplot,weightplot)
                plt.xlabel('Energy (keV)')
                plt.title("Sampled source spectrum"+ str(sumW))
                plt.show()
            
            if self.spectrumFromXls:
                rightUnitScale=1
                if self.source_dict["energyUnit"]=="eV":
                    rightUnitScale=0.001
                elif self.source_dict["energyUnit"]=="MeV":
                    rightUnitScale=1000
                
                spectrum=[]
                colEnergy=None
                colFluence=None
                i=0
                for sh in xlrd.open_workbook(self.source_dict["pathXlsSpectrum"]).sheets():
                    foundData=0
                    for row in range(sh.nrows):
                        for col in range(sh.ncols):
                            myCell = sh.cell(row, col)
        #                    print("\n\n Sample materials : ",self.myMaterials[imat])
                            if myCell.value == self.source_dict["energyColumnKey"]:
                                colEnergy=col
                                startRow=row
                                foundData+=1
                            if myCell.value == self.source_dict["fluenceColumnKey"]:
                                colFluence=col
                                foundData+=1
                            
                        if foundData==2:
                            break
                        
                    if colEnergy==None:
                        raise Exception (f'Energy column key {self.source_dict["energyColumnKey"]} not found in the xls file')
                    if colFluence==None:
                        raise Exception (f'Energy column key {self.source_dict["fluenceColumnKey"]} not found in the xls file')
                        
                        
                    sumFluence=0
                    for row in range(startRow+1, sh.nrows):
                        myCellFluence = sh.cell(row, colFluence)
                        fluence=myCellFluence.value
                        myCellEnergy = sh.cell(row, colEnergy)
                        energy=myCellEnergy.value* rightUnitScale
                        if myCell.value is not None:
                            spectrum.append([energy,fluence])
                            sumFluence+=fluence
                            
                arraySpectrum=np.asarray(spectrum)
                
                # plt.figure()
                # plt.plot(arraySpectrum[:,0],arraySpectrum[:,1])
                # plt.xlabel('Energy (keV)')
                # plt.title("Source spectrum")
                # plt.show()
                
                #re-sampling at sourcre energy sampling (to get faster calculation at the end)
                energyplot=[]
                weightplot=[]
                den=spectrum[1][0]-spectrum[0][0]
                    
                Nen=len(spectrum)
                Nbin=int((spectrum[-1][0]-spectrum[0][0])//self.source_dict["myEnergySampling"])
                # Nbin=int(np.ceil(Nen/self.source_dict["myEnergySampling"]/2))
                n=0
                totWeight=0
                for i in range(Nbin-1):
                    currBin=0
                    weightBin=0
                    energyBin=0
                    while currBin<self.source_dict["myEnergySampling"]:
                        weightBin+=spectrum[n][1]
                        energyBin+=spectrum[n][1]*spectrum[n][0]
                        n+=1
                        currBin=currBin+den
                    if weightBin!=0:
                        # self.mySpectrum.append((energyBin/weightBin,weightBin))
                        energyplot.append(energyBin/weightBin)
                        weightplot.append(weightBin)
                    totWeight+=weightBin
                    
                currBin=0
                weightBin=0
                energyBin=0
                while n<Nen:
                    weightBin+=spectrum[n][1]
                    energyBin+=spectrum[n][1]*spectrum[n][0]
                    n+=1
                if weightBin!=0:
                    # self.mySpectrum.append((energyBin/weightBin,weightBin))
                    energyplot.append(energyBin/weightBin)
                    weightplot.append(weightBin)
                
                totFlux=0
                for i in range(len(energyplot)):
                    flux=weightplot[i]/totWeight
                    totFlux+=flux
                    if flux>0.001:
                        self.mySpectrum.append((energyplot[i],flux))
                
                
                finalSpectrum=self.mySpectrum
                # k=0
                # while self.mySpectrum[0][1]/totWeight<0.001:
                #     self.mySpectrum.pop(0)
                #     energyplot.pop(0)
                #     weightplot.pop(0)
                #     k+=1
                
                plt.figure()
                plt.plot(energyplot,weightplot)
                plt.xlabel('Energy (keV)')
                plt.title(f"Xls source spectrum total flux {totFlux}")
                plt.show()
                    
        return
            
            
        raise ValueError("type of source not recognized")
        
    
    def getText(self,node):
        return node.childNodes[0].nodeValue
    
    
if __name__ == "__main__":
    
    # name = 'Water'
    # comment = 'Defined by chemical formula'
    # composition='H2O'
    # sp.Spek.make_matl(matl_name=name, matl_density=1.0, chemical_formula=composition,matl_comment=comment)
    
    
    s = sp.Spek(kvp=20,th=12, targ='Mo', dk=0.5)
    s.filter('Air', 800)
    # s.filter('Water', 20)
    # s.filter('C', 0.2)
    # s.filter('Cu', 0.05)
    spectrum=s.get_spectrum(flu=True)
    
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
    print("Flux", sumFlux*1e-8)
    print("Energy stdDev:", stdDev)
    
    plt.figure()
    plt.plot(spectrum[0],spectrum[1]/np.max(spectrum[1]))
    plt.xlabel('Energy (keV)')
    # plt.savefig('/Users/quenot/Library/Mobile Documents/com~apple~CloudDocs/Thèse/ComparaisonSimu/W40kVp.png')
    plt.show()
    
    
    
    
    
