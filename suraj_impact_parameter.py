# python suraj_impact_parameter.py 'qso catalogue' 'galaxy catalogue'
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 13:52:15 2018

@author: surajpoudel
"""


import sys
from math import *
from astropy.io import ascii
from astropy.table import Table
import numpy as np
import math
from math import cos, sin, radians
if len(sys.argv) == 4:
    qso_cat=sys.argv[1]
    gal_cat=sys.argv[2]
    b=int(sys.argv[3])
data1=Table.read(qso_cat, format='ascii') #read the quasar files
data2=Table.read(gal_cat, format='ascii') #read the galaxy files with name, ra, dec , and redshift
#data1=Table.read('qso.txt', format='ascii') #read the quasar files
#data2=Table.read('galaxy.txt', format='ascii') #read the galaxy files with name, ra, dec , and redshift
H0=70 # Hubble's constant
WM=0.3 # Matter density
WV=0.7 # Vacuum density
row1=data1['col1']
name1=data1['col2']
iden1=data1['col3']
Ra1=data1['col4'] 
Dec1=data1['col5']
type1=data1['col6']
zq=data1['col8']
mag=data1['col9']
name2=data2['col1']
Ra2=data2['col2']
Dec2=data2['col3']
z=data2['col4']


#final=[]
for r,n1,id1,t1,z1,mg,i,k in zip(row1,name1,iden1,type1,zq,mag,Ra1,Dec1):
     #print n1
     ff=[r,n1,id1,t1,z1,i,k]
     for n2,j,l,z2 in zip(name2,Ra2,Dec2,z):
             #if (cos(radians(90-k))*cos(radians(90-l))+sin(radians(90-k))*sin(radians(90-l))*cos(radians(i-j)))< cos(radians(8/60)):
             #m = ((((i-j)*cos(radians((k+l)/2)))**2+(k-l)**2)**(1/2))
             s = (((i-j)*cos((radians)((k+l)/2)))**2+(k-l)**2)
             #s = (((i-j)*cos(radians(l/2)))**2+(k-l)**2)
             sep=math.sqrt(s)          # seperation in degrees
             # initialize constants
                 
             WR = 0.        # Omega(radiation)
             WK = 0.        # Omega curvaturve = 1-Omega(total)
             c = 299792.458 # velocity of light in km/sec
             Tyr = 977.8    # coefficent for converting 1/H into Gyr
             DTT = 0.5      # time from z to now in units of 1/H0
             DTT_Gyr = 0.0  # value of DTT in Gyr
             age = 0.5      # age of Universe in units of 1/H0
             age_Gyr = 0.0  # value of age in Gyr
             zage = 0.1     # age of Universe at redshift z in units of 1/H0
             zage_Gyr = 0.0 # value of zage in Gyr
             DCMR = 0.0     # comoving radial distance in units of c/H0
             DCMR_Mpc = 0.0
             DCMR_Gyr = 0.0
             DA = 0.0       # angular size distance
             DA_Mpc = 0.0
             DA_Gyr = 0.0
             kpc_DA = 0.0
             DL = 0.0       # luminosity distance
             DL_Mpc = 0.0
             DL_Gyr = 0.0   # DL in units of billions of light years
             V_Gpc = 0.0
             a = 1.0        # 1/(1+z), the scale factor of the Universe
             az = 0.5       # 1/(1+z(object))
             h = H0/100.
             WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528
             WK = 1-WM-WR-WV
             az = 1.0/(1+1.0*z2)
             age = 0.
             n=1000         # number of points in integrals
             for i1 in range(n):
                a = az*(i1+0.5)/n
                adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
                age = age + 1./adot
             zage = az*age/n
             zage_Gyr = (Tyr/H0)*zage
             DTT = 0.0
             DCMR = 0.0

             # do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
             for i2 in range(n):
                 a = az+(1-az)*(i2+0.5)/n
                 adot = sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
                 DTT = DTT + 1./adot
                 DCMR = DCMR + 1./(a*adot)
        
             DTT = (1.-az)*DTT/n
             DCMR = (1.-az)*DCMR/n
             age = DTT+zage
             age_Gyr = age*(Tyr/H0)
             DTT_Gyr = (Tyr/H0)*DTT
             DCMR_Gyr = (Tyr/H0)*DCMR
             DCMR_Mpc = (c/H0)*DCMR

             # tangential comoving distance

             ratio = 1.00
             x = sqrt(abs(WK))*DCMR
             if x > 0.1:
                if WK > 0:
                   ratio =  0.5*(exp(x)-exp(-x))/x
                else:
                   ratio = sin(x)/x
             else:
                   y = x*x
                   if WK < 0: y = -y
                   ratio = 1. + y/6. + y*y/120.
             DCMT = ratio*DCMR
             DA = az*DCMT
             DA_Mpc = (c/H0)*DA
             kpc_DA = DA_Mpc/206.264806
             DA_Gyr = (Tyr/H0)*DA
             DL = DA/(az*az)
             DL_Mpc = (c/H0)*DL
             DL_Gyr = (Tyr/H0)*DL

             sep_as=sep*3600          #separation in arc sec
             impact=sep_as*kpc_DA     # impact parameter in kpc
             #print z
             if impact<b: # if the impact parameter is less than b.0 kpc
                 
                  print r,n1,id1,t1,z1,i,k,n2,j,l,z2,sep,impact # print the output 
                  #print r,n1,id1,t1,z1,i,k,n2,j,l,z2,sep,impact
   # read the output as: qso_row, qso_name, identity, type, qso_z, ra_qso, dec_qso, galaxy_name, 
   # ra_gal, dec_gal, galaxy_z, seperation in degrees, impact parameter in kpc
            
