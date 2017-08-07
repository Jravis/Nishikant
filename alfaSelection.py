"""
In this code we are creating catalog for alfa-alfa ra and dec for sdss crossid
and than we try to create random catalog using that
"""
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ndarray
import math as m
from scipy.integrate import quad

def E(x, Om, OL):
    return 1.0 / (Om * ((1.0 + x) ** 3.0) + OL) ** 0.5


Om = 0.3
OL = 0.7

dH = (30000.0 / 7.0)  # Hubble distance


def DL(red):
    dm, err = quad(E, 0, red, args=(Om, OL))
    dl = (dm * (1 + red) * dH)
    return dl


path = os.path.abspath('/home/mavrick/Documents/data_stacking/AlFA70/AlfaSDSSDR7_sandy_7690_0.csv')
path1 = os.path.abspath('/home/mavrick/Documents/data_stacking/AlFA70/AlfasdssDR7/alfa40SDD.txt')
path2 = os.path.abspath('/home/mavrick/Documents/data_stacking/AlFA70/AlfasdssDR7/crossid71_sandy_7690.csv')


objid = np.genfromtxt(path, usecols=(0), delimiter=',', dtype=None)
ra = np.genfromtxt(path, usecols=(1), delimiter=',', dtype=None)
dec = np.genfromtxt(path, usecols=(2), delimiter=',', dtype=None)
g = np.genfromtxt(path, usecols=(4), delimiter=',', dtype=None)
r = np.genfromtxt(path, usecols=(5), delimiter=',', dtype=None)
absMI = np.genfromtxt(path, usecols=(16), delimiter=',', dtype=None)
absMg = np.genfromtxt(path, usecols=(17), delimiter=',', dtype=None)
Z = np.genfromtxt(path, usecols=(8), delimiter=',', dtype=None)

#objid1 = np.genfromtxt(path2, usecols=(4),delimiter=',', dtype=None)
objid1 = np.genfromtxt(path1, usecols=(3), dtype=None)





print len(Z), len(ra), len(objid1)
Ra = []
Dec = []
G = []
R = []
MI = []
Mg = []
redshift = []

"""
for i in range(len(objid)):
    if objid[i] in objid1:
        Ra.append(ra[i])
        Dec.append(dec[i])
        G.append(g[i])
        R.append(r[i])
        MI.append(absMI[i])
        Mg.append(absMg[i])
        redshift.append(Z[i])

print len(redshift), len(Ra)

plt.figure(1)
plt.scatter(Ra, Dec, c=redshift)
plt.colorbar()
plt.show()
"""