#import numpy as np
#import scipy as sp
from scipy import constants as conts

def ret_energy(wavelength):
    return (conts.h * conts.c) / (wavelength*conts.e) / conts.angstrom

def ret_wl(energy):
    return (conts.h * conts.c) / (energy*conts.e) / conts.angstrom