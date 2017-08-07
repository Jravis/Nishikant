from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ndarray
import math as m
from scipy.integrate import quad
from mpl_toolkits.mplot3d import Axes3D
def E(x, Om, OL):
    return 1.0 / (Om * ((1.0 + x) ** 3.0) + OL) ** 0.5


Om = 0.3
OL = 0.7

dH = (30000.0 / 7.0)  # Hubble distance


def DL(red):
    dm, err = quad(E, 0, red, args=(Om, OL))
    dl = (dm * (1 + red) * dH)
    return dl


path = os.path.abspath('/home/mavrick/Documents/data_stacking/AlFA70/AlfasdssDR7/test.txt')
path1 = os.path.abspath('/home/mavrick/Documents/data_stacking/AlFA70/AlfasdssDR7/alfa40SDD.txt')
path2 = os.path.abspath('/home/mavrick/Documents/data_stacking/AlFA70/AlfasdssDR7/AlfaOCproperties_sandy_7690_0.csv')


objid = np.genfromtxt(path2, usecols=0, delimiter=',', dtype=None)

u = np.genfromtxt(path2, usecols=3, delimiter=',', dtype=None)
g = np.genfromtxt(path2, usecols=4, delimiter=',', dtype=None)
absMu = np.genfromtxt(path2, usecols=18, delimiter=',', dtype=None)




"""
objid1 = np.genfromtxt(path1, usecols=(3), dtype=None)
AGCNr = np.genfromtxt(path1, usecols=0, dtype=None)



agcNr = np.genfromtxt(path, usecols=0,delimiter='\t', dtype=None)
mHI = np.genfromtxt(path, usecols=6, delimiter='\t', dtype=None)
Vhelio = np.genfromtxt(path, usecols=8, delimiter='\t', dtype=None)
W50 = np.genfromtxt(path, usecols=9, delimiter='\t', dtype=None)
errW50 = np.genfromtxt(path, usecols=10, delimiter='\t', dtype=None)

data = ndarray((len(objid), 4), float)

for i in range(len(objid)):
    for j in range(len(objid1)):
        if objid[i] == objid1[j]:
            data[i, 0] = AGCNr[j]
            data[i, 1] = absMu[i]
            data[i, 2] = g[i]
            data[i, 3] = u[i]

Data = ndarray((len(objid), 7), float)

for j in range(data[:, 0].size):
    for i in range(len(agcNr)):
        if agcNr[i] == data[j, 0]:
            Data[j, 0] = mHI[i]
            Data[j, 1] = W50[i]
            Data[j, 2] = data[j, 1]
            Data[j, 3] = Vhelio[i]
            Data[j, 4] = errW50[i]
            Data[j, 5] = data[j, 2]
            Data[j, 6] = data[j, 3]

np.savetxt('test1.csv', zip(Data[:, 0], Data[:, 1], Data[:, 2], Data[:, 3], Data[:, 4], Data[:, 5], Data[:, 6]), delimiter=',', fmt='%f', header='M_HI, W50, absMu, Vhelio, errW50, g, u')

"""

MHI = np.genfromtxt('test1.csv', usecols=0, delimiter=',', dtype=None)
w50 = np.genfromtxt('test1.csv', usecols=1, delimiter=',', dtype=None)
Mu = np.genfromtxt('test1.csv', usecols=2, delimiter=',', dtype=None)
Vhelio = np.genfromtxt('test1.csv', usecols=3, delimiter=',', dtype=None)
errW50 = np.genfromtxt('test1.csv', usecols=4, delimiter=',', dtype=None)

g = np.genfromtxt('test1.csv', usecols=5, delimiter=',', dtype=None)
u = np.genfromtxt('test1.csv', usecols=6, delimiter=',', dtype=None)

print len(errW50)

a = []
b = []
c = []
d = []
U = []
G = []
for i in range(len(MHI)):
    if Vhelio[i] >= 1000:
        if MHI[i] != 0.000000:
            a.append(MHI[i])
        if w50[i] != 0.000000:
            b.append(w50[i])
            d.append(errW50[i])
        if Mu[i] != 0.000000:
            if Mu[i] != -9999.00000:
                c.append(Mu[i])
            else:
                c.append(0.0)
        if u[i] != 0.000000:
            U.append(u[i])
        if g[i] != 0.000000:
            G.append(g[i]-u[i])


print len(a), len(b), len(c), len(d), len(G), len(G)
plt.figure(1)



fig = plt.figure(1)
ax = fig.add_subplot(111,projection='3d')
ax.scatter(c, G, a, zdir=u'z', s=10, c='blue')
plt.ylim(-5.0, 2.0)


"""
plt.hist2d(c, a, bins=100, norm=LogNorm())
plt.xlabel("absM_u")
plt.ylabel("log(HI)")
plt.colorbar()
plt.figure(2)
plt.scatter(a, b)
plt.ylabel("W50")
plt.xlabel("log(HI)")
plt.yscale("log")

plt.figure(3)
plt.hist2d(a, b, bins= 100, norm=LogNorm())
plt.ylabel("W50")
plt.xlabel("log(HI)")
plt.colorbar()
plt.figure(4)
plt.scatter(c, a)
plt.xlabel("absMu")
plt.ylabel("log(HI)")



plt.figure(5)
plt.scatter(c, b)
plt.errorbar(c, b, yerr=d, linestyle="None")
plt.xlabel("absM_u")
plt.ylabel("W50")
plt.yscale("log")
"""
plt.show()

