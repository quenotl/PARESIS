#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 15:57:17 2021

@author: quenot
"""
import numpy as np
import xlrd 

     
def getDeltaBeta(self, sourceSpectrum):
    energyRange=[5,100]
#        print("Materials :", self.myMaterials)
    pathTablesDeltaBeta ='Samples/DeltaBeta/TablesDeltaBeta.xls'
    deltaBetaDoc=xlrd.open_workbook(pathTablesDeltaBeta)
    
    i=0
    for sh in xlrd.open_workbook(pathTablesDeltaBeta).sheets():
        for imat in range(len(self.myMaterials)):
            delta=[]
            beta=[]
            if self.myMaterials[imat]=="None":
                for energy in range(energyRange[0],energyRange[1]):
                    delta.append((energy,0.))
                    beta.append((energy,0.))
                self.beta.append(beta)
                self.delta.append(delta)
                pass
                
            for col in range(sh.ncols):
                row=0
                myCell = sh.cell(row, col)
#                    print("\n\n Sample materials : ",self.myMaterials[imat])
                if myCell.value == self.myMaterials[imat]:
                    
                    row=row+3
                    for energy in range(energyRange[0],energyRange[1]):
                        j=i
                        while j==i:
                            currentCell=sh.cell(row+1,col)
                            if float(currentCell.value)>=energy*1e3:
                                delta.append((energy,float(sh.cell(row,col+1).value)))
                                beta.append((energy,float(sh.cell(row,col+2).value)))
                                i=i+1
                                
                            row=row+1
                        
                    self.beta.append(beta)
                    self.delta.append(delta)
                    
#                        print("delta :", delta)
        a, b, c=np.shape(self.delta)
        if a<len(self.myMaterials):
            raise ValueError("One or more materials have not been found in delta beta tables")
        return