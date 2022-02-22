#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 09:16:10 2020

@author: quenot
"""

import numpy as np
from scipy.constants import c, e, h


def getk(energy):
    """Get the wavenumber from energy in eV

    Parameters
    ----------
    energy : array_like
        _description_

    Returns
    -------
    _type_
        _description_
    """
    k = 2*np.pi*energy*e/(h*c)
    return k


if __name__ == "__main__":
    K = getk(25000)
    print("k=", K)
