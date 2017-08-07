import numpy as np 
import matplotlib.pyplot as plt
import scipy
import random 
from mpl_toolkits.mplot3d import Axes3D
from numpy import ndarray
import math as m
from scipy.integrate import quad


Om = 0.3
OL = 0.7
cm_to_Mpc = 3.0857e24

def E(x, Om, OL):
    return 1.0/(Om*((1.0+x)**3.0)+OL)**(0.5)

dH = (30000.0/7.0)
dm, err = quad(E, 0.0, 0.05, args=(Om,OL))
region = (dm*1.05*dH)

dm, err = quad(E, 0.0, 0.0023, args=(Om, OL))
region2 = (dm*1.0023*dH)
print region, region2

#random.seed(8388608)
#random.seed(87889091)
random.seed(81100289)
count1 = 0

"""

#------------------------------------------------------
#Conversion from HH:MM:SS to degree

RA = 0.0
DEC = 0.0
i = "Yes"
while (i !="No"):
    HH = float(raw_input("Enter HH for RA :- "))
    MM = float(raw_input("Enter MM for RA :- "))
    SS = float(raw_input("Enter SS for RA :- "))
    
    dd = float(raw_input("Enter dd for DEC :-"))
    M1M1 = float(raw_input("Enter MM for DEC :-"))
    S1S1 = float(raw_input("Enter SS for DEC :-"))
    
    RA = (HH*15.0)+(MM/4.0)+(SS/240.0)
    if (dd < 0.0):
        dd = -1.0*dd
        DEC = -1.0*(dd + (M1M1*5.0/3.0 + S1S1*5.0/180.0)/100.0)
    else:    
        DEC = (dd + (M1M1*5.0/3.0 + S1S1*5.0/180.0)/100.0)

    print RA , DEC

    print "You want to continue if no type No"
    i = raw_input("Enter:- ")
#-----------------------------------------------------
"""


x1 = ndarray((12190,), float)
y1 = ndarray((12190,), float)
z1 = ndarray((12190,), float)
ra = ndarray((12190,), float)
dec = ndarray((12190,), float)
dist = ndarray((12190,), float)

while (count1< 12190):

    theta = np.arccos(2*np.random.uniform(0, 1)-1)
    phi = np.random.uniform(0, 2*np.pi)
    u = np.random.uniform(0, 1)
    r = u**(1.0/3.0)
    if r*region >= region2:
        if ((112.5*m.pi/180.0) <= phi):
            if(phi <= (247.5*m.pi/180)):
                if ((0.0*m.pi/180) <= ((90*np.pi/180)-theta)<= (m.pi*18.0/180)):
                    ra[count1] = phi*180/m.pi
                    dec[count1]= (90-theta*180/m.pi)
                    z1[count1] = r*region*m.cos(theta)
                    y1[count1] = r*region*m.sin(theta)*m.sin(phi)
                    x1[count1] = r*region*m.sin(theta)*m.cos(phi)
                    dist[count1] = r*region
                    count1+=1
                if((23.9*np.pi/180)<= ((90*np.pi/180)-theta)<= (30.2*np.pi/180)):
                        ra[count1] = phi*180/m.pi
                        dec[count1]= (90-theta*180/m.pi)
                        z1[count1] = r*region*m.cos(theta)
                        y1[count1] = r*region*m.sin(theta)*m.sin(phi)
                        x1[count1] = r*region*m.sin(theta)*m.cos(phi)
                        dist[count1] = r*region
                        count1+=1

        if ((129.5*m.pi/180.0) <= phi):
            if(phi <= (232.0*m.pi/180)):
                if ((18.0*m.pi/180) <= ((90*np.pi/180)-theta) <= (m.pi*22.4/180)):
                        ra[count1] = phi*180/m.pi
                        dec[count1]= (90-theta*180/m.pi)
                        z1[count1] = r*region*m.cos(theta)
                        y1[count1] = r*region*m.sin(theta)*m.sin(phi)
                        x1[count1] = r*region*m.sin(theta)*m.cos(phi)
                        dist[count1] = r*region
                        count1+=1


#------------------------------------------------------------------------------------------------
"""

    if ((330*m.pi/180.0) <= phi):
        if(phi <= (360.0*m.pi/180)):
            if ((0.0*m.pi/180) <= theta ):# or
                if((theta<= (m.pi*2.4/180))):# or
                    ra[count1] = phi*180/m.pi
                    dec[count1]= theta*180/m.pi
                    z1[count1] = r*region*m.cos((90.0*m.pi/180)-theta)
                    y1[count1] = r*region*m.sin((90.0*m.pi/180)-theta)*m.sin(phi)
                    x1[count1] = r*region*m.sin((90.0*m.pi/180)-theta)*m.cos(phi)
                    count1+=1 
            if ((6.0*m.pi/180) <= theta ):# or
                if((theta<= (m.pi*10/180))):# or
                    ra[count1] = phi*180/m.pi
                    dec[count1]= theta*180/m.pi
                    z1[count1] = r*region*m.cos((90.0*m.pi/180)-theta)
                    y1[count1] = r*region*m.sin((90.0*m.pi/180)-theta)*m.sin(phi)
                    x1[count1] = r*region*m.sin((90.0*m.pi/180)-theta)*m.cos(phi)
                    count1+=1 

            if((13.8*np.pi/180)<= theta):
                if(theta <= (16.0*np.pi/180)):
                    ra[count1] = phi*180/m.pi
                    dec[count1]= theta*180/m.pi
                    z1[count1] = r*region*m.cos((90.0*m.pi/180)-theta)
                    y1[count1] = r*region*m.sin((90.0*m.pi/180)-theta)*m.sin(phi)
                    x1[count1] = r*region*m.sin((90.0*m.pi/180)-theta)*m.cos(phi)
                    count1+=1 
            if((21*np.pi/180)<= theta):
                if(theta <= (32.0*np.pi/180)):
                    ra[count1] = phi*180/m.pi
                    dec[count1]= theta*180/m.pi
                    z1[count1] = r*region*m.cos((90.0*m.pi/180)-theta)
                    y1[count1] = r*region*m.sin((90.0*m.pi/180)-theta)*m.sin(phi)
                    x1[count1] = r*region*m.sin((90.0*m.pi/180)-theta)*m.cos(phi)
                    count1+=1 
"""
#--------------------------------------------------------------------------------------
"""
    if ((0.0*m.pi/180.0) <= phi):
        if(phi <= (45.0*m.pi/180)):
            if ((0.0*m.pi/180) <= theta ):# or
                if((theta<= (m.pi*2.4/180))):# or
                    ra[count1] = phi*180/m.pi
                    dec[count1]= theta*180/m.pi
                    z1[count1] = r*region*m.cos((90.0*m.pi/180)-theta)
                    y1[count1] = r*region*m.sin((90.0*m.pi/180)-theta)*m.sin(phi)
                    x1[count1] = r*region*m.sin((90.0*m.pi/180)-theta)*m.cos(phi)
                    count1+=1 
            if ((6.0*m.pi/180) <= theta ):# or
                if((theta<= (m.pi*10/180))):# or
                    ra[count1] = phi*180/m.pi
                    dec[count1]= theta*180/m.pi
                    z1[count1] = r*region*m.cos((90.0*m.pi/180)-theta)
                    y1[count1] = r*region*m.sin((90.0*m.pi/180)-theta)*m.sin(phi)
                    x1[count1] = r*region*m.sin((90.0*m.pi/180)-theta)*m.cos(phi)
                    count1+=1 

            if((13.8*np.pi/180)<= theta):
                if(theta <= (16.0*np.pi/180)):
                    ra[count1] = phi*180/m.pi
                    dec[count1]= theta*180/m.pi
                    z1[count1] = r*region*m.cos((90.0*m.pi/180)-theta)
                    y1[count1] = r*region*m.sin((90.0*m.pi/180)-theta)*m.sin(phi)
                    x1[count1] = r*region*m.sin((90.0*m.pi/180)-theta)*m.cos(phi)
                    count1+=1 
            if((21*np.pi/180)<= theta):
                if(theta <= (32.0*np.pi/180)):
                    ra[count1] = phi*180/m.pi
                    dec[count1]= theta*180/m.pi
                    z1[count1] = r*region*m.cos((90.0*m.pi/180)-theta)
                    y1[count1] = r*region*m.sin((90.0*m.pi/180)-theta)*m.sin(phi)
                    x1[count1] = r*region*m.sin((90.0*m.pi/180)-theta)*m.cos(phi)
                    count1 += 1
    

"""

np.savetxt("ranCatCorrd1.txt", zip(x1, y1, z1, ra, dec), fmt='%f', delimiter='\t', header='x,y,z,ra,dec',newline='\n')


Ra_HI = np.genfromtxt('test.txt', usecols=(1), delimiter='\t', dtype=None)
Dec_HI = np.genfromtxt('test.txt', usecols=(2), delimiter='\t', dtype=None)

print count1
"""
fig = plt.figure(1)
ax = fig.add_subplot(111,projection ='3d')
ax.scatter(x1,y1,z1, zdir=u'z', s = 10, c='blue')
"""
fig = plt.figure(2)
plt.plot(x1, y1, "bo", ms= 1)
fig = plt.figure(3)
plt.plot(x1, z1, "bo", ms= 1)

fig = plt.figure(4)
plt.plot(y1, z1, "bo", ms= 1)
fig = plt.figure(5)
plt.plot(ra, dec, "bo", ms= 3)
plt.plot(Ra_HI, Dec_HI, "go", ms= 3)

plt.figure(6)
n_bins = 50
n, bins, patches = plt.hist(dist, n_bins, histtype='step')
plt.xlabel("Distance(Mpc)")
plt.ylabel("N")
plt.title("Variation of Number of galaxies with Distance")

plt.show()
