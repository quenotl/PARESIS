#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 09:16:10 2020

@author: quenot
"""

import numpy as np


def getk(energy):
    """
    energy in eV
    """
    h=6.62607015e-34
    c=2.99792458e8
    e=1.60217663e-19
    k=2*np.pi*energy*e/(h*c)
    return k



     
if __name__ == "__main__":
    K=getk(25000)
    print("k=", K)
