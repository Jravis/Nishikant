

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.ascii as ascii
import matplotlib.gridspec as gridspec

data1 = ascii.read('/home/sandeep/Nishi/data/a40.datafile1.csv', delimiter=',', guess=False)

agcNr = data1['AGCNr']
ra_HI = data1['RAdeg_HI']
dec_HI = data1['Decdeg_HI']
M_HI = data1['logMsun']
D = data1['Dist']
Vhelio = data1['Vhelio']
W50 = data1['W50']
HICode = data1['HIcode']
HIflux = data1['HIflux']
HIflux = np.asarray(HIflux)

agcNr = np.asarray(agcNr)  # Alfa id
ra_HI = np.asarray(ra_HI)  # HI ra
dec_HI = np.asarray(dec_HI)  # HI dec
M_HI = np.asarray(M_HI)   #  HI Galaxy mass in log
Dist = np.asarray(D)  # HI galaxy distance in Mpc
Vhelio = np.asarray(Vhelio)  # HI galaxy heliocentric velocity
W50 = np.asarray(W50)  # HI line width W50 parameter
HICode = np.asarray(HICode)

Vlimit = 1000.0  # to remove
index = (Vhelio > Vlimit) * (HICode == 1)
print index
print HICode[index]
print Vhelio[index]

c = 3e5
z = Vhelio[index]/c
"""
plt.figure(1, figsize=(15,6))
gs = gridspec.GridSpec(1, 2)
ax1 = plt.subplot(gs[0, 0])
ax1.plot(z, M_HI[index], "ro", ms = 3, label=r"$V_{helio}$ > 1000, HICODE = 1")

plt.xlabel(r'$Z$', fontsize='large', fontstyle='italic', weight='extra bold')
plt.ylabel(r'$M_{HI}$', fontsize='large', fontstyle='italic', weight='extra bold')
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=14)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=14)
ax1.legend(loc=4)
ax2 = plt.subplot(gs[0, 1])
ax2.plot(Dist[index], M_HI[index], "bo", ms = 3, label=r"$V_{helio}$ > 1000, HICODE = 1")
plt.xlabel(r'$D$', fontsize='large', fontstyle='italic', weight='extra bold')
plt.ylabel(r'$log_{10}(M_{HI})$', fontsize='large', fontstyle='italic', weight='extra bold')
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=14)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=14)
ax2.legend(loc=4)

plt.tight_layout()
"""

from decimal import *

N_jk, xedges, yedges = np.histogram2d(M_HI[index], np.log10(W50[index]), bins=8)
print N_jk.shape
N_jk = N_jk.T

"""
plt.figure(7)  # , figsize=(10, 5))
twod_hist = plt.hist2d((M_HI[index]), np.log10(W50[index]), bins=40, cmap='RdBu')  # , norm=LogNorm())
plt.colorbar()
plt.figure(8)  # , figsize=(10, 5))
plt.imshow(N_jk, interpolation='nearest', origin='low', cmap='RdBu')
plt.colorbar()
"""
print np.diff(xedges)

delta_M = np.diff(xedges)[0]
delta_W = np.diff(yedges)[0]
print delta_W

bin_arr_M = []
bin_arr_W50 = []
for i in xrange(len(xedges)):
    if i + 1 < len(xedges):
        tmp = xedges[i]
        tmp1 = xedges[i + 1]
        arry = [float(Decimal('%0.6f' % tmp)), float(Decimal('%0.6f' % tmp1))]
        bin_arr_M.append(arry)

        tmp = yedges[i]
        tmp1 = yedges[i + 1]
        arry = [float(Decimal('%0.6f' % tmp)), float(Decimal('%0.6f' % tmp1))]
        bin_arr_W50.append(arry)

bin_M = np.asarray(bin_arr_M)
bin_W50 = np.asarray(bin_arr_W50)

print bin_M.shape
print ""
print bin_W50.shape
log_W = np.log10(W50[index])

print ""
from scipy.integrate import quad, dblquad

HI_mass = M_HI[index]

H_ijk = np.zeros((len(HIflux[index]), len(xedges), len(yedges)), dtype=np.float64)

HI_flux = HIflux[index]
print log_W
S50 = 0.0

for i in xrange(len(HI_flux)):
    for j in xrange(len(bin_M)):
        for k in xrange(len(bin_W50)):

            if np.min(bin_M[j]) <= HI_mass[i] <= np.max(bin_M[j]):
                if np.min(bin_W50[k]) <= log_W[i] <= np.max(bin_W50[k]):

                    if log_W[i] < 2.5:
                        S50 = 0.5 * log_W[i] - 1.14 - 0.130
                    if np.log10(HI_flux[i]) > S50:
                        C = 1.
                    else:
                        C = 0.

                    if log_W[i] >= 2.5:
                        S50 = log_W[i] - 2.39 - 0.130
                    if np.log10(HI_flux[i]) > S50:
                        C = 1.
                    else:
                        C = 0.

                    H_ijk[i, j, k] = C

a = len(HI_flux)
phi = np.zeros((len(xedges), len(yedges)), dtype=np.float64)
phi1 = np.ones((len(xedges), len(yedges)), dtype=np.float64)
norm = np.zeros(a, dtype=np.float64)


for nn in xrange(0, 1):
    for j in xrange(len(bin_M)):
        for k in xrange(len(bin_W50)):
            Sum1 = 0.0
            for i in xrange(len(HI_flux)):
                Sum = 0.0
                for l in xrange(len(bin_M)):
                    for m in xrange(len(bin_W50)):
                        Sum += H_ijk[i, l, m]*phi1[l, m]
                        #print H_ijk[i, l, m], phi1[l, m]
               # norm[i] = Sum
            #for i in xrange(len(HIflux[index])):
#                Sum1 += H_ijk[i, j, k]/Sum
                print ""
                print Sum, H_ijk[i, j, k]
            #print N_jk[j, k], Sum1
           # phi[j, k] = N_jk[j, k]/Sum1

            #phi1 = phi

#print phi




