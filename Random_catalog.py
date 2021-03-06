
import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.cosmology import FlatLambdaCDM
import astropy.io.ascii as ascii
from astropy.constants import c  # the speed of light
from astroML.plotting import hist
import scipy
import TwoPCF
from astroML.correlation import two_point


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
    if type(ra) != float and  type(ra)!=np.float64:
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

Count=0
for i in xrange(len(ra_HI)):
    if 135 < ra_HI[i] < 230:
        if 0 < dec_HI[i] < 18:
            Count+=1
print ""
print "Emanuell dataset"
print Count


index = (M_HI > 9.9)*(W50 > 18) * (HICode == 1)
#index = (HICode == 1)

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

print "Emanuell dataset"
print len(W50_50)
#===============================================================

mass = ['10.5-10.0', '10.0-9.5', '9.5-8.5']
fname = '/dataspace/sandeep/ALFALFA_data/selectionFunc_HImass_bins/slecFunc_%s_bin.txt' % mass[1]

Dist_psi = np.genfromtxt(fname, usecols=0)
psi_slection = np.genfromtxt(fname, usecols=1)

fint_psi = scipy.interpolate.interp1d(Dist_psi, psi_slection, kind='cubic')


fname = '/dataspace/sandeep/ALFALFA_data/rfi_selection/rfi_completeness.txt'
RFI_helio = np.genfromtxt(fname, usecols=0, delimiter=',')
RFI_selec = np.genfromtxt(fname, usecols=1, delimiter=',')


fint_RFI = scipy.interpolate.interp1d(RFI_helio, RFI_selec, kind='cubic')


index1 = (M_HI_50 <= 10.0) * (M_HI_50 >= 9.5)
mass = M_HI_50[index1]

data_ra = ra_HI_50[index1]
data_dec = dec_HI_50[index1]
data_z = z_cmb_50[index1]
data_dist = Dist_50[index1]
np.savetxt('/dataspace/sandeep/Nishi_plots/Nishi_dataCat.txt', zip(data_ra, data_dec, data_z),
           fmt='%0.6e\t%0.6e\t%0.6e', delimiter='\t')
data = radec_to_xyz(data_ra, data_dec)
data_D = np.zeros((len(data_ra), 3), dtype=np.float64)
data_D[:, 0] = np.multiply(data[0, :], data_dist)
data_D[:, 1] = np.multiply(data[1, :],data_dist)
data_D[:, 2] = np.multiply(data[2, :], data_dist)


#===== Random Catalog====

ncount = 0
rand_selection = []
count1 = 0
count2 = 0
count3 = 0
rand_dist = []
rand_ra = []
rand_dec = []

np.random.seed(7337777)

while ncount < 5*len(mass):
    key = 'None'
    #while key == 'None':
    u = np.random.uniform(0, 1)
    r = u**(1.0/3.0)
        #temp_dist = np.random.uniform(min(Dist_psi), max(Dist_psi))
        #selection = fint_psi(temp_dist)
    temp_dist = r*max(Dist_psi)
        #if temp_dist >= min(Dist_psi):
        #    selection = fint_psi(temp_dist)
        #    u1 = np.random.uniform(min(psi_slection), 1.0)
         #   if u1 <= selection:
         #       key = 'Found'

#    theta = (np.random.uniform(0, 1))*15.0+72.0
    #Ra = np.random.uniform(0, 1)*95+135

    phi = np.random.uniform(0, 2*np.pi)
    theta = np.arccos(2*np.random.uniform(0, 1)-1)
    
    Dec = np.subtract(90, np.degrees(theta))
    Ra = np.degrees(phi)
    if 135 < Ra < 230:
        if 0.0 < Dec < 18.0:

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
rand_ra = np.asarray(rand_ra)
rand_dec = np.asarray(rand_dec)

V_cmb = Ho * rand_dist
Z_cmb = V_cmb / c_km_s

np.savetxt('/dataspace/sandeep/Nishi_plots/Nishi_randomCat.txt', zip(rand_ra, rand_dec, Z_cmb),
           fmt='%0.6e\t%0.6e\t%0.6e', delimiter='\t')
data = radec_to_xyz(rand_ra, rand_dec)
data_R = np.zeros((len(rand_ra), 3), dtype=np.float64)
data_R[:, 0] = np.multiply(data[0, :], rand_dist)
data_R[:, 1] = np.multiply(data[1, :], rand_dist)
data_R[:, 2] = np.multiply(data[2, :], rand_dist)



#==== 2PCF Part ====
#bins = np.linspace(0.0, 21.0, 50)
bins = np.logspace(np.log10(0.5), np.log10(20.0), 15)
xi = TwoPCF.two_point(data_D, bins, method='standard', data_R=data_R, random_state=None)

# PLotting Part      

"""
plt.figure(1, figsize=(8, 6))
h1 = hist(v_cmb_50, bins='knuth', histtype='stepfilled', ec='k', fc='#AAAAAA')
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=14)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=14)
plt.ylabel('N', fontsize='large', fontstyle='italic', weight='extra bold')
plt.xlabel(r'$v_{cmb}$', fontsize='xx-large', fontstyle='italic', weight='extra bold')
plt.figure(2, figsize=(8, 6))
h1 = hist(z_cmb_50, bins='knuth', histtype='stepfilled', ec='k', fc='#AAAAAA')
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=14)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=14)
plt.ylabel('N', fontsize='large', fontstyle='italic', weight='extra bold')
plt.xlabel(r'$z_{cmb}$', fontsize='xx-large', fontstyle='italic', weight='extra bold')

bin1 = 10**np.linspace(np.log10(min(rand_dist)), np.log10(max(rand_dist)), 500)
hist1, xedges = np.histogram(rand_dist, bins=bin1)

plt.figure(3, figsize=(10, 4))
h2 = hist(Dist_cmb_50, bins=50, histtype='stepfilled', ec='k', fc='#AAAAAA', label='data')
#h3 = hist(rand_dist, bins=bin1, histtype='step', ec='g', fc='#AAAAAA', label='random')
plt.plot(xedges[:-1], hist1, 'r-', linewidth=1,label='random' )
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=14)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=14)
plt.ylabel('N', fontsize='large', fontstyle='italic', weight='extra bold')
plt.xlabel(r'$D_{cmb}$', fontsize='xx-large', fontstyle='italic', weight='extra bold')
plt.legend()
plt.xlim(10, 220)

plt.savefig('/dataspace/sandeep/Nishi_plots/redshift_dist.eps', dpi=100)


plt.figure(4, figsize=(9, 9))
plt.plot(ra_HI_50, dec_HI_50, 'g.', ms =2)
plt.plot(rand_ra, rand_dec, 'r.', ms=2)
plt.minorticks_on()
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=14)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=14)
plt.xlabel(r'$Ra_{HI}$', fontsize='large', fontstyle='italic', weight='extra bold')
plt.ylabel(r'$Dec_{HI}$', fontsize='large', fontstyle='italic', weight='extra bold')
plt.savefig('/dataspace/sandeep/Nishi_plots/rand_data.eps', dpi=100)


plt.figure(6)
plt.plot(data_R[:, 0], data_R[:, 1], 'bo', ms=2)
plt.tick_params(axis='both', which='minor', length=5, width=2, labelsize=14)
plt.tick_params(axis='both', which='major', length=8, width=2, labelsize=14)
plt.xlabel('x', fontsize='large', fontstyle='italic', weight='extra bold')
plt.ylabel('y', fontsize='large', fontstyle='italic', weight='extra bold')

plt.savefig('/dataspace/sandeep/Nishi_plots/xy.eps', dpi=100)

plt.figure(7)
plt.plot(data_R[:, 0], data_R[:, 2], 'go', ms=2)
plt.xlabel('x', fontsize='large', fontstyle='italic', weight='extra bold')
plt.ylabel('z', fontsize='large', fontstyle='italic', weight='extra bold')

plt.savefig('/dataspace/sandeep/Nishi_plots/xz.eps', dpi=100)


plt.figure(8)
plt.plot(data_R[:, 1], data_R[:, 2], 'ro', ms=2)
plt.xlabel('y', fontsize='large', fontstyle='italic', weight='extra bold')
plt.ylabel('z', fontsize='large', fontstyle='italic', weight='extra bold')
plt.savefig('/dataspace/sandeep/Nishi_plots/yz.eps', dpi=100)
plt.figure(9)
plt.plot(rand_dist, rand_selection, 'ro', ms=2)
plt.xlabel(r'$Dist_{random}$', fontsize='large', fontstyle='italic', weight='extra bold')
plt.ylabel(r'$f(x)$', fontsize='large', fontstyle='italic', weight='extra bold')

plt.savefig('/dataspace/sandeep/Nishi_plots/rand_selection.eps', dpi=100)
"""

x = bins[:-1]
plt.figure(10)
plt.loglog(bins[:-1], xi, 'r-o', ms=5, linewidth=2)
plt.loglog(x, x**-1.1, 'b-', linewidth=2)


plt.show()


