"""
 In particular, we select galaxies over a contiguous
rectangular sky region of 1 700 deg2 (135◦ < RA < 230◦ and 0◦ < Dec < 18◦≈and in the redshift range z ≈ 0.0023 − 0.05
(vCMB = 700 − 15 000 km s−1 ).
The sample is restricted to “Code 1” ALFALFA detections, i.e. it is comprised only by confidently detected
extragalactic sources (S/NHI > 6.5). In addition, parent sample sources have a
combination of observed 21cm flux (SHI ) and 21cm lineprofile width (W50) that
places them in the region of the {SHI , W50}–plane where the completeness of the
ALFALFA survey is at least 50% [see Sec. 6 and Fig. 12 in Haynes et al., 2011].
Lastly,
 the sample is limited to linewidths W50 > 18 km s−1 and HI masses2
M HI > 107.5 M⊙ .

"""

import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.cosmology import FlatLambdaCDM
import astropy.io.ascii as ascii
from astropy.constants import c  # the speed of light
from astroML.plotting import hist
import scipy

cosmo = FlatLambdaCDM(70.0, 0.3, Neff=0)
Ho = 70.0  #
c_km_s = c.to('km/s').value

def hms_to_deg(h, m, s):
    return 15. * (h + m / 60. + s / 3600.)


def sdms_to_deg(sign, d, m, s):
    sign = 1. - 2. * (sign == '-')
    return sign * (d + m / 60. + s / 3600.)


def jd_to_mjd(t):
    return t - 2400000.5


def mag_to_flux(m, me, zp):
    """Convert magnitude and magnitude error to flux, given a zeropoint."""

    f = 10. ** (0.4 * (zp - m))
    fe = math.log(10.) * 0.4 * me * f

    return f, fe


def radec_to_xyz(ra, dec):
    if type(ra) != float:
        temp = np.zeros((3, len(ra)), dtype=np.float64)

        temp[0, :] = np.sin(np.deg2rad(dec)) * np.cos(np.deg2rad(ra))
        temp[1, :] = np.sin(np.deg2rad(dec)) * np.sin(np.deg2rad(ra))
        temp[2, :] = np.cos(np.deg2rad(dec))
        return temp
    else:
        xx = np.sin(np.deg2rad(dec)) * np.cos(np.deg2rad(ra))
        yy = np.sin(np.deg2rad(dec)) * np.sin(np.deg2rad(ra))
        zz = np.cos(np.deg2rad(dec))
        return np.array([xx, yy, zz], dtype=np.float64)


def cmb_dz(ra, dec):
    """See http://arxiv.org/pdf/astro-ph/9609034
     CMBcoordsRA = 167.98750000 # J2000 Lineweaver
     CMBcoordsDEC = -7.22000000
    """

    # J2000 coords from NED
    CMB_DZ = 371000. / 299792458.
    CMB_RA = 168.01190437
    CMB_DEC = -6.98296811
    CMB_XYZ = radec_to_xyz(CMB_RA, CMB_DEC)

    coords_xyz = radec_to_xyz(ra, dec)

    dz = CMB_DZ * np.dot(CMB_XYZ, coords_xyz)

    return dz


def helio_to_cmb(z, ra, dec):
    """Convert from heliocentric redshift to CMB-frame redshift.

    Parameters
    ----------
    z : float
        Heliocentric redshift.
    ra, dec: float
        RA and Declination in degrees (J2000).
    """

    dz = -cmb_dz(ra, dec)
    one_plus_z_pec = np.sqrt((1. + dz) / (1. - dz))
    one_plus_z_CMB = (1. + z) / one_plus_z_pec

    return one_plus_z_CMB - 1.


def cmb_to_helio(z, ra, dec):
    """Convert from CMB-frame redshift to heliocentric redshift.

    Parameters
    ----------
    z : float
        CMB-frame redshift.
    ra, dec: float
        RA and Declination in degrees (J2000).
    """

    dz = -cmb_dz(ra, dec)
    one_plus_z_pec = math.sqrt((1. + dz) / (1. - dz))
    one_plus_z_helio = (1. + z) * one_plus_z_pec

    return one_plus_z_helio - 1.


"""
def random_catalog(Nsample, selection_func, rfi_func, dist_slection, rfi_helio):
    :param Nsample:
    :param selection_func:
    :param rfi_func:
    :param dist_slection:
    :param rfi_helio:
    :return:

#    ra_rand = np.zeros(20*Nsample, dtype=np.float64)
#    dec_rand = np.zeros(20*Nsample, dtype=np.float64)
#    dist_rand = np.zeros(20*Nsample, dtype=np.float64)
#    selection_rand = np.zeros(20*Nsample, dtype=np.float64)

    ra_rand = []
    dec_rand = []
    dist_rand = []
    selection_rand = []

    ncount = 0


    for i in xrange(0, 50):
        u = random.uniform(0, 1.0001)
        if min(selection_func) <= u <= max(selection_func):
            j = 0
            while u <= psi_slection[j]:
                j += 1

            zi_cmb = min(dist_slection)*cosmo.H0.value/c_km_s
            #Vol1 = cosmo.comoving_volume(dist_slection[j-1]*cosmo.H0.value/c_km_s).value - cosmo.comoving_volume(zi_cmb).value
            Vol1 = coVol(dist_slection[j-1]*cosmo.H0.value/c_km_s) - coVol(zi_cmb)
            ngal = int(Vol1*Ngal)
            print ngal


            while ncount < ngal:
                tht = np.arccos(2 * random.uniform(0, 1) - 1)
                phi = random.uniform(0, 2 * np.pi)
                r = random.uniform(0, 1.0001)
                if min(dist_slection) <= r*dist_slection[j-1] <= dist_slection[j-1]:
                    if 135 < np.degrees(phi) < 230:
                        if 3.5 < 90-np.degrees(tht) < 16.5:
                            V_cmb = Ho*r*dist_slection[j-1]
                            Z_cmb = V_cmb/c_km_s
                            Z_helio = cmb_to_helio(Z_cmb, m.degrees(phi), 90.-m.degrees(tht))
                            V_helio = Z_helio*c_km_s
                            k = 0
                            while V_helio < rfi_helio[k]:
                                k += 1

                            if rfi_func[k-1] >= 0.8:

                                dist_rand.append(r * dist_slection[j - 1])
                                ra_rand.append(np.degrees(phi))
                                dec_rand.append((90 - np.degrees(tht)))
                                selection_rand.append(psi_slection[j])
                                ncount += 1



    ra_rand = np.asarray(ra_rand)
    dec_rand = np.asarray(dec_rand)
    dist_rand = np.asarray(dist_rand)
    selection_rand = np.asarray(selection_rand)

    xyz = radec_to_xyz(ra_rand, dec_rand)
    xyz[0, :] = np.multiply(xyz[0, :], dist_rand)
    xyz[1, :] = np.multiply(xyz[1, :], dist_rand)
    xyz[2, :] = np.multiply(xyz[2, :], dist_rand)

    return ra_rand, dec_rand, xyz, dist_rand

    """
#  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

data1 = ascii.read('/home/sandeep/Nishi/data/a40.datafile1.csv', delimiter=',', guess=False)
ra_HI = data1['RAdeg_HI']
dec_HI = data1['Decdeg_HI']
M_HI = data1['logMsun']
D = data1['Dist']
Vhelio = data1['Vhelio']
W50 = data1['W50']
HICode = data1['HIcode']
HIflux = data1['HIflux']
RMS = data1['RMS']
SNR = data1['SNR']


ra_HI = np.asarray(ra_HI)  # HI ra
dec_HI = np.asarray(dec_HI)  # HI dec
M_HI = np.asarray(M_HI)  #  HI Galaxy mass in log
Dist = np.asarray(D)  # HI galaxy distance in Mpc
Vhelio = np.asarray(Vhelio)  # HI galaxy heliocentric velocity
W50 = np.asarray(W50)  # HI line width W50 parameter
HICode = np.asarray(HICode)
HIflux = np.asarray(HIflux)
RMS = np.asarray(RMS)
SNR = np.asarray(SNR)

index = (M_HI > 7.5)*(W50 > 18) * (HICode == 1)

ra_HI = ra_HI[index]  # HI ra
dec_HI = dec_HI[index]  # HI dec
M_HI = M_HI[index]  #  HI Galaxy mass in log
Dist = Dist[index]  # HI galaxy distance in Mpc
Vhelio = Vhelio[index]  # HI galaxy heliocentric velocity
W50 = W50[index]  # HI line width W50 parameter
HICode = HICode[index]
HIflux = HIflux[index]

z_helio = Vhelio/c_km_s
z_cmb = helio_to_cmb(z_helio, ra_HI, dec_HI)
v_cmb = z_cmb*c_km_s
Dist_cmb = v_cmb/Ho

red_index = (v_cmb >= 700.)*(v_cmb <= 15000)

ra_HI = ra_HI[red_index]  # HI ra
dec_HI = dec_HI[red_index]  # HI dec
M_HI = M_HI[red_index]  #  HI Galaxy mass in log
Dist = Dist[red_index]  # HI galaxy distance in Mpc
Vhelio = Vhelio[red_index]  # HI galaxy heliocentric velocity
W50 = W50[red_index]  # HI line width W50 parameter
HICode = HICode[red_index]
HIflux = HIflux[red_index]
v_cmb = v_cmb[red_index]
Dist_cmb = Dist_cmb[red_index]
z_cmb = z_cmb[red_index]

# above 50% completness of ALFALFA

log_W = np.log10(W50)

ra_HI_50 = []
dec_HI_50 = []
M_HI_50 = []
Dist_50 = []
Vhelio_50 = []
W50_50 = []
HICode_50 = []
HIflux_50 = []
v_cmb_50 = []
Dist_cmb_50 = []
z_cmb_50 = []

for i in xrange(len(log_W)):
    if log_W[i] < 2.5:
        S50_1 = 0.5 * log_W[i] - 1.14 - 0.130
        if np.log10(HIflux[i]) >= S50_1:
            if 135 < ra_HI[i] < 230:
                if 0 < dec_HI[i] < 18:
                    ra_HI_50.append(ra_HI[i])
                    dec_HI_50.append(dec_HI[i])
                    M_HI_50.append(M_HI[i])
                    Dist_50.append(Dist[i])
                    Vhelio_50.append(Vhelio[i])
                    W50_50.append(W50[i])
                    HICode_50.append(HICode[i])
                    HIflux_50.append(HIflux[i])
                    v_cmb_50.append(v_cmb[i])
                    Dist_cmb_50.append(Dist_cmb[i])
                    z_cmb_50.append(z_cmb[i])

    if log_W[i] >= 2.5:
        S50_2 = log_W[i] - 2.39 - 0.130
        if np.log10(HIflux[i]) >= S50_2:
            if 135 < ra_HI[i] < 230:
                if 0 < dec_HI[i] < 18:
                    ra_HI_50.append(ra_HI[i])
                    dec_HI_50.append(dec_HI[i])
                    M_HI_50.append(M_HI[i])
                    Dist_50.append(Dist[i])
                    Vhelio_50.append(Vhelio[i])
                    W50_50.append(W50[i])
                    HICode_50.append(HICode[i])
                    HIflux_50.append(HIflux[i])
                    v_cmb_50.append(v_cmb[i])
                    Dist_cmb_50.append(Dist_cmb[i])
                    z_cmb_50.append(z_cmb[i])


print "test"
ra_HI_50 = np.asarray(ra_HI_50, dtype=np.float64)
dec_HI_50 = np.asarray(dec_HI_50, dtype=np.float64)
M_HI_50 = np.asarray(M_HI_50, dtype=np.float64)
Dist_50 = np.asarray(Dist_50, dtype=np.float64)
Vhelio_50 = np.asarray(Vhelio_50, dtype=np.float64)
W50_50 = np.asarray(W50_50, dtype=np.float64)
HICode_50 = np.asarray(HICode_50, dtype=np.float64)
HIflux_50 = np.asarray(HIflux, dtype=np.float64)
v_cmb_50 = np.asarray(v_cmb_50, dtype=np.float64)
Dist_cmb_50 = np.asarray(Dist_cmb_50, dtype=np.float64)
z_cmb_50 = np.asarray(z_cmb_50, dtype=np.float64)


mass = ['10.5-10.0', '10.0-9.5', '9.5-8.5']
fname = '/dataspace/sandeep/ALFALFA_data/selectionFunc_HImass_bins/slecFunc_%s_bin.txt' % mass[1]
#fname = '/dataspace/sandeep/ALFALFA_data/selectionFunc_mass_thres/test_3.txt'

Dist_psi = np.genfromtxt(fname, usecols=0)
psi_slection = np.genfromtxt(fname, usecols=1)

fint_psi = scipy.interpolate.interp1d(Dist_psi, psi_slection, kind='cubic')


fname = '/dataspace/sandeep/ALFALFA_data/rfi_selection/rfi_completeness.txt'
RFI_helio = np.genfromtxt(fname, usecols=0, delimiter=',')
RFI_selec = np.genfromtxt(fname, usecols=1, delimiter=',')


fint_RFI = scipy.interpolate.interp1d(RFI_helio, RFI_selec, kind='cubic')


index1 = (M_HI_50 <= 10.0) * (M_HI_50 >= 9.5)
mass = M_HI_50[index1]


#theta = (np.random.uniform(0, 1, size=5*len(mass)))*15.0+72.0
#rand_dec = np.subtract(90, theta)
#rand_ra = np.random.uniform(0, 1, size=5*len(mass))*95+135

#np.random.seed(7715177)
#temp_dist = np.random.uniform(min(Dist_psi), max(Dist_psi), size=5*len(mass))
#selection = fint(temp_dist)

ncount = 0
rand_selection = []
count1 = 0
count2 = 0
count3 = 0
rand_dist = []
rand_ra = []
rand_dec = []
rand_dist = []
print Dist_psi

np.random.seed(7337777)
while ncount < 20*len(mass):
    key = 'None'
    while key == 'None':
        temp_dist = np.random.uniform(min(Dist_psi), max(Dist_psi))
        selection = fint_psi(temp_dist)
        u = np.random.uniform(min(psi_slection), 1.0)
        if u <= selection:
            key = 'Found'

    theta = (np.random.uniform(0, 1))*15.0+72.0
    Dec = np.subtract(90, theta)
    Ra = np.random.uniform(0, 1)*95+135

    V_cmb = Ho*temp_dist

    Z_cmb = V_cmb/c_km_s
    Z_helio = cmb_to_helio(Z_cmb, Ra, Dec)
    V_helio = Z_helio*c_km_s

    rfi_selection = fint_RFI(V_helio)
    rfi = np.random.uniform(0, 1.0)
    if rfi <= rfi_selection:
        rand_ra.append(Ra)
        rand_dec.append(Dec)
        rand_dist.append(temp_dist)
        rand_selection.append(u)
        ncount += 1

rand_dist = np.asarray(rand_dist)
data = radec_to_xyz(rand_ra, rand_dec)
data[0, :] = np.multiply(data[0, :], rand_dist)
data[1, :] = np.multiply(data[1, :], rand_dist)
data[2, :] = np.multiply(data[2, :], rand_dist)


plt.figure(1, figsize=(8, 6))
h1 = hist(v_cmb_50, bins='knuth', histtype='stepfilled', ec='k', fc='#AAAAAA')
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=14)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=14)
plt.ylabel(r'$N$', fontsize='large', fontstyle='italic', weight='extra bold')
plt.xlabel(r'$v_{cmb}$', fontsize='xx-large', fontstyle='italic', weight='extra bold')

plt.figure(2, figsize=(8, 6))
h1 = hist(z_cmb_50, bins='knuth', histtype='stepfilled', ec='k', fc='#AAAAAA')
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=14)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=14)
plt.ylabel(r'$N$', fontsize='large', fontstyle='italic', weight='extra bold')
plt.xlabel(r'$z_{cmb}$', fontsize='xx-large', fontstyle='italic', weight='extra bold')

bin1 = 10**np.linspace(np.log10(min(rand_dist)), np.log10(max(rand_dist)), 500)
hist1, xedges = np.histogram(rand_dist, bins=bin1)

plt.figure(3, figsize=(10, 4))
h2 = hist(Dist_cmb_50, bins=50, histtype='stepfilled', ec='k', fc='#AAAAAA', label='data')
#h3 = hist(rand_dist, bins=bin1, histtype='step', ec='g', fc='#AAAAAA', label='random')
plt.plot(xedges[:-1], hist1, 'g-', linewidth=1)
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=14)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=14)
plt.ylabel(r'$N$', fontsize='large', fontstyle='italic', weight='extra bold')
plt.xlabel(r'$D_{cmb}$', fontsize='xx-large', fontstyle='italic', weight='extra bold')
plt.legend()
plt.xlim(10, 220)
plt.figure(4, figsize=(8, 6))
plt.plot(ra_HI_50, dec_HI_50, 'go', ms =3)
plt.plot(rand_ra, rand_dec, 'ro', ms=3)
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=14)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=14)
plt.xlabel(r'$Ra_{HI}$', fontsize='large', fontstyle='italic', weight='extra bold')
plt.ylabel(r'$Dec_{HI}$', fontsize='large', fontstyle='italic', weight='extra bold')


plt.figure(6)
plt.plot(data[0, :], data[1, :], 'bo', ms=2)
plt.figure(7)
plt.plot(data[0, :], data[2, :], 'go', ms=2)
plt.figure(8)
plt.plot(data[1, :], data[2, :], 'ro', ms=2)

plt.figure(9)
plt.plot(rand_dist, rand_selection, 'ro', ms=2)

plt.show()
