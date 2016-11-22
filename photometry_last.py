'''
    
    Photometry SUPER LOTIS Analysis
    by Zeynep Yasemin Kalender
    
    
    '''


import os, csv
from glob import glob
from astropy.io import fits #For opening fits files
from astropy.coordinates import SkyCoord #For creating sky coordinates
from astropy import units as u #For assigning units to numbers
from photutils import aperture_photometry, SkyCircularAperture ,SkyCircularAnnulus #For perfomring photometry
from astropy.cosmology import WMAP9 as cosmo
from astropy.table import Table
from astropy.io import ascii
from astropy.table import hstack
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time


'''Type the name of the supernova when you are asked, also be careful to change the coordinate
    In this code , the coordinate below is for SN13DAD 
    In order to check the code, one can try with sn13dad '''

sn_name = raw_input('enter the supernova name: >')



#Path Names

super_lotis_path = '/Users/zeynepyaseminkalender/Documents/MyFiles_Pyhton/SuperLOTIS_final'
date_search_path = os.path.join(super_lotis_path, '13*') #this path choses all 13* directories.



def photometry():
    
    filters = ['I','B','V','R'] #filter names
   
   # create a table here for Photometry results
    photo_table=Table(names=('Residuals','Error','Time'))
    Photometry = open('photometry_'+str(sn_name)+'.txt','w')
    
    for f in filters:
        
        #create a txt file where the Photometry Results table will be saved
        
        
        #join the two paths given above
        search_str = os.path.join(date_search_path,str(sn_name)+str(f)+'.fits')
    
        #names will be the supernova fits files for different filters
        for name in glob(search_str):
        
            with fits.open(name) as analysis:
                
                N = analysis[0].header['NCOMBINE'] # number of images stacked in the fit file
                final_date = convert_time(analysis) # The Time series must be MJD
                
                #Position of Supernova as coordinate, aperture is calculated with the coordinate
                cordinate = SkyCoord('01:48:08.66 +37:33:29.22', unit = (u.hourangle, u.deg))
                #coordinates = SkyCoord('01:48:08.66 +37:33:29.22',unit = (u.hourangle, u.deg))
                
                #exp_time= analysis[0].header['EXPTIME']
                #data_error = np.sqrt(N*analysis[0].data)
                
                aperture = SkyCircularAperture(cordinate, r=5*u.arcsec)
                annulus_apertures = SkyCircularAnnulus(cordinate, r_in=6*u.arcsec, r_out=8*u.arcsec)
                
                aperture_area = np.pi * 3 ** 2
                annulus_area = np.pi * (8 ** 2 - 6 ** 2)
                
                rawflux_table = aperture_photometry(analysis[0], aperture)
                bkgflux_table = aperture_photometry(analysis[0], annulus_apertures )
                
                residual_sum = rawflux_table[0]['aperture_sum'] - bkgflux_table[0]['aperture_sum'] * aperture_area / annulus_area
                error=np.sqrt(np.abs(residual_sum)*N)
                photo_table.add_row((residual_sum,error,final_date))
                
    print "***** Photometry Results without Calibration ****"
    print photo_table
    print >> Photometry, photo_table
    Photometry.close()
    calibration(N)

def star_photometry():
    
    filters = ['I','B','V','R'] #filter names
    
    # create a table here for Photometry results
    photo_table=Table(names=('Residuals','Error','Time'))
    S_Photometry = open('star_'+str(sn_name)+'.txt','w')
    
    for f in filters:
        
        
        #join the two paths given above
        search_str = os.path.join(date_search_path,str(sn_name)+str(f)+'.fits')
        
        #names will be the supernova fits files for different filters
        for name in glob(search_str):
            
            with fits.open(name) as analysis:
                
                N = analysis[0].header['NCOMBINE'] # This is for calculation how many images stacked in the fit file
                final_date = convert_time(analysis) # The Time series must be MJD
                
                '''Position of Supernova as coordinate, aperture is calculated with the coordinate '''
                cordinate = SkyCoord('01:47:29.318 +37:26:32.12', unit = (u.hourangle, u.deg))
                coordinates = SkyCoord('01:47:29.318 +37:26:32.12',unit = (u.hourangle, u.deg))
                
                #exp_time= analysis[0].header['EXPTIME']
                #data_error = np.sqrt(N*analysis[0].data)
                
                aperture = SkyCircularAperture(cordinate, r=5*u.arcsec)
                annulus_apertures = SkyCircularAnnulus(coordinates, r_in=6*u.arcsec, r_out=8*u.arcsec)
                
                aperture_area = np.pi * 3 ** 2
                annulus_area = np.pi * (8 ** 2 - 6 ** 2)
                
                rawflux_table = aperture_photometry(analysis[0], aperture)
                bkgflux_table = aperture_photometry(analysis[0], annulus_apertures )
                
                residual_sum = rawflux_table[0]['aperture_sum'] - bkgflux_table[0]['aperture_sum'] * aperture_area / annulus_area
                error=np.sqrt(np.abs(residual_sum)*N)
                photo_table.add_row((residual_sum,error,final_date))



    print >> S_Photometry, photo_table
    S_Photometry.close()



def convert_time(analysis):
    
    #analysis is the fits file in the photometry function
    t = Time(analysis[0].header['DATE-OBS'])
    t.format = 'mjd'
    result = t.value
    return result





def calibration(number):
    
    #Opening our photometry results that we calculated,
    #....and star results I just calculated in Some other code
    
    cal_sn = ascii.read('photometry_'+str(sn_name)+'.txt','r')
    cal_star = ascii.read('star_'+str(sn_name)+'.txt','r')
    results=Table()
    mag_star=13.018
    
    results["Count_SN"] = cal_sn["Residuals"]
    results["Count_Star"] = cal_star["Residuals"]
    results["Inst_SN"] = -2.5*np.log10(cal_sn["Residuals"])
    results["Inst_Star"] = -2.5*np.log10(cal_star["Residuals"])
    results["Z_point"] = mag_star-results["Inst_Star"]
    results["Mag_Star"] = results["Inst_Star"]+results["Z_point"]
    results["Mag_SuperN"] = -2.5*np.log10(cal_sn["Residuals"])+results["Z_point"]
    
    #file["Error_mag"] = 1.086*number/file['Mag_sn_'+str(f)]
    print " ***** Calibrated Results ****"
    print results
    
        
    '''
        plt.errorbar(file["Time"],file['Mag_sn_'+str(f)],file["Error_mag"],file["Error_mag"]*0,'o',label='SN')
        plt.plot(file["Time"],file["Mag_Star"],'k*',label='STAR')
        plt.xlabel("MJD_time")
        plt.ylabel("Aparent_Mag")
        plt.title("Sn13dad_"+str(f))
        plt.legend()
        plt.savefig("Mag_"+str(sn_name)+str(f))
        #plt.show()
        print >> results, file
    
    results.close()
    '''


#calling the functions
star_photometry()
photometry()















