# This code for analysing lUMINOSITY DISTANCE-REDSHIFT RELATIONSHIP
# using astropy packages.

# Date= 10 JUNE 2016
# Author = Z. Yasemin Kalender



#table  and plot packages
from astropy.io import ascii
import numpy as  np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits
# cosmology packages
from astropy.cosmology import WMAP9 as cosmo
from astropy.cosmology import  FlatwCDM , FlatLambdaCDM
import astropy.units as u
from astropy.cosmology import Planck13, z_at_value

##### Some Cosmological Constants..
H0=cosmo.H(0)
cosmo = FlatwCDM(name='SNLS3+WMAP7', H0=71.58, Om0=0.262, w0=-1.016)


# FRIEDMAN 2015 PAPER TABLE 9

path= '/Users/zeynepyaseminkalender/Documents/MyFiles_Pyhton/week5'

#obtaining Data needed

Friedman_data9= Table.read('friedman_table9.fit')

print(Friedman_data9)


#redshift value
z_values=Friedman_data9["zCMB"]

#Luminosity distance calculation
lum_dis=cosmo.luminosity_distance(z_values)
lum_dis=lum_dis.value
print " Luminosity-distance after calculated:" , lum_dis


#Luminosity distance module calculation

lum_dist_modulus=cosmo.distmod(z_values)
print "Luminosity distance modulus values are :" ,lum_dist_modulus
lum_dist_modulus=lum_dist_modulus.value
print "Luminosity distance modulus values are :" ,lum_dist_modulus

# Method 1 without eliminating Mpc



plt.figure(1)
plt.errorbar(z_values,lum_dis,0*Friedman_data9["e_zCMB"],Friedman_data9["e_zCMB"],fmt='o')
#plt.scatter(z_values,lum_dis)
plt.xlabel("Red Shift")
plt.ylabel("Luminosity Distance - Mpc ")
plt.title("Friedman 2015 Paper-Luminosity distance versus Z ")
plt.show()
plt.close()

plt.figure(2)
plt.errorbar(z_values,lum_dist_modulus,0*Friedman_data9["e_zCMB"],Friedman_data9["e_zCMB"],fmt='o')
plt.xlabel("Red Shift")
plt.ylabel("Luminosity Distance Modulus - Magnitude")
plt.title("Friedman 2015 Paper-Luminosity distance modulus versus Z ")
plt.show()



'''
# Method 2  to eliminate Mpc, converting into string
#converting to Mpc we can multimpy by #tmp=tmp*3.09*10**22
tmp=str(lum_dis)
tmp=tmp.replace("Mpc","")
print (" Luminosity Distance: ")
print tmp
print (" Red shift: ")
print z_values
#Plot
plt.scatter(z_values,tmp)
plt.xlabel("Red Shift")
plt.ylabel("Luminosity Distance - Mpc ")
plt.Title("Friedman 2015 Paper-Luminosity distance versus Z ")
plt.show()
'''


# Checking the number of variable

print "Luminosity distance Row:",len(lum_dis)
print " Red Shift Row:" , len(z_values)
#print " New Variable tmp lenth: " , len(tmp)
#print "type of tmp:" , type(tmp)
print "type of z:" ,type(z_values)











