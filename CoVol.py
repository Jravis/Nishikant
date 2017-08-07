"Routine for Comoving volume"
import numpy as np
from scipy.integrate import quad
from astropy.cosmology import FlatLambdaCDM
import random
from astropy.constants import c  # the speed of light
from astroML.plotting import hist
import matplotlib.pyplot as plt
random.seed(49390927)
cosmo = FlatLambdaCDM(70.0, 0.3, Neff=0)

phi = 0.0044
alpha = -1.35
Mstar = 10**9.98
shecterFun = lambda t: np.log(10)*phi*(t/Mstar)**(alpha+1)*np.exp(-(t/Mstar))/t

Ngal, err = quad(shecterFun, 10.0**9.5, 10.0**10) #, args=())
print Ngal


def SolidAngle(red):
    #return np.pi*(1.22*0.21*(1+red)*0.5/45)**2 # GMRT
    return np.pi*(0.21*(1+red)*0.5/305)**2  # ALFALFA

# ALFALFA survey region 1700 deg**2


def DV1(z2):
    Da = cosmo.angular_diameter_distance(z2).value
    Ez = cosmo.H(z2).value/cosmo.H0.value
    dH = c/cosmo.H0.value
    return (SolidAngle(z2)*dH*(1+z2)**2*Da**2)/Ez


def DV2(z):
    return cosmo.differential_comoving_volume(z).value * 1700*np.deg2rad(1.)**2  #SolidAngle(z)


def coVol(red):
    z = red
#    dm1, err1 = quad(DV1, 0.0, z)
    dm2, err2 = quad(DV2, 0.0, z)
    #return dm1, dm2
    return dm2


mass = ['10.5-10.0', '10.0-9.5', '9.5-8.5']
fname = '/dataspace/sandeep/ALFALFA_data/selectionFunc_HImass_bins/slecFunc_%s_bin.txt' % mass[1]
#fname = '/dataspace/sandeep/ALFALFA_data/selectionFunc_mass_thres/test_3.txt'

print fname
dist_slection = np.genfromtxt(fname, usecols=0)
selection_func = np.genfromtxt(fname, usecols=1)
c_km_s = c.to('km/s').value

fname = '/dataspace/sandeep/ALFALFA_data/rfi_selection/rfi_completeness.txt'
print fname
RFI_helio = np.genfromtxt(fname, usecols=0, delimiter=',')
RFI_selec = np.genfromtxt(fname, usecols=1, delimiter=',')

count=0
dist = []
#for i in xrange(0, 50):


#    zi_cmb = min(dist_slection)*cosmo.H0.value/c_km_s
#    Vol1 = coVol(dist_slection[j-1]*cosmo.H0.value/c_km_s) - coVol(zi_cmb)
#    ngal = int(Vol1*Ngal)
#    print min(dist_slection), (dist_slection[j-1])
#    print ngal
"""
while count <= 10000:
    u = random.uniform(min(selection_func), max(selection_func))
    j = 0
    if 0.999914 < u <= max(selection_func):
        j = int(random.uniform(0, 64))
    else:
        while u <= selection_func[j]:
            j += 1

    print j, u
    r = random.uniform(0, 1.0)
    if min(dist_slection) <= r * dist_slection[j - 1] <= dist_slection[j - 1]:
        dist.append(r*dist_slection[j-1])
        count += 1

"""


for i in xrange(0, 50):
    zi_cmb = min(dist_slection)*cosmo.H0.value/c_km_s

    j = int(random.uniform(0, len(dist_slection)))
    Vol1 = coVol(dist_slection[j]*cosmo.H0.value/c_km_s) - coVol(zi_cmb)
    ngal = int(Vol1*Ngal)

    print j, ngal, dist_slection[j]
    while count <= ngal:
        r = random.uniform(0, 1.0)
        if min(dist_slection) <= r * dist_slection[j] <= dist_slection[j]:
            dist.append(r*dist_slection[j])
            count += 1



plt.figure(3, figsize=(10, 6))
#h1 = hist(Dist_cmb_50, bins='knuth', histtype='step', ec='k', fc='#AAAAAA', label='Knuth')# 26 bins
h3 = hist(dist, bins=100, histtype='step', ec='g', fc='#AAAAAA', label='Const')

plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=14)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=14)
plt.ylabel(r'$N$', fontsize='large', fontstyle='italic', weight='extra bold')
plt.xlabel(r'$D_{cmb}$', fontsize='xx-large', fontstyle='italic', weight='extra bold')
plt.legend()
plt.yscale('log')
plt.show()