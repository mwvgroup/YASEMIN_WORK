
from astropy.io import fits #For opening fits files
from astropy.coordinates import SkyCoord #For creating sky coordinates
from astropy import units as u #For assigning units to numbers
from photutils import aperture_photometry, SkyCircularAperture , CircularAperture#For perfomring photometry
from astropy.cosmology import WMAP9 as cosmo
from astropy.table import Table
from astropy.io import ascii
from glob import glob
import os.path ,os ,csv
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time


'''Type the name of the supernova when you are asked, also be careful to change the coordinate
    In this code , the coordinate below is for SN13DAD 
    In order to check the code, one can try with sn13dad '''
sn_name = raw_input('enter the supernova name: >')

'''Path names'''
super_lotis_path = '/Users/zeynepyaseminkalender/Documents/MyFiles_Pyhton/SuperLOTIS_final'
date_search_path = os.path.join(super_lotis_path, '13*') #this path choses all 13* directories.


'''Here is the Photometry code starts'''
def photometry():
    
    filters = ['I','B','V','R'] #filter names
    
    for f in filters:
        
        '''create a table where the results will be attached'''
        All_data = Table(names=('Count','Error','Time'))
        
        '''create a txt file where the Photometry Results will be saved'''
        file = open('photometry_'+str(sn_name)+str(f)+'.txt','w')
        
        search_str = os.path.join(date_search_path,str(sn_name)+str(f)+'.fits') #join the two pathes given above
    
        for name in glob(search_str): #names will be the supernova fits files for different filters
        
            with fits.open(name) as analysis:
                
                N = analysis[0].header['NCOMBINE'] # This is for calculation how many images stacked in the fit file
                final_date = convert_time(analysis) # The Time series must be MJD
                
                '''Position of Supernova as coordinate, aperture is calculated with the coordinate '''
                cordinate = SkyCoord('01:48:08.66 +37:33:29.22', unit = (u.hourangle, u.deg))
                aperture = SkyCircularAperture(cordinate, r=3*u.arcsec)
                
                '''Error Calculation '''
                exp_time= analysis[0].header['EXPTIME'] # error calculation
                data_error = np.sqrt(N*analysis[0].data)
                
                '''Performing Photometry and recording them in a --tbl-- '''
                tbl=Table(aperture_photometry(analysis[0],aperture,error=data_error),names=('Count','Count_Error','x_center','y_center','center_input'))
                #tbl.keep_columns(['Count','Count_Error'])
                
                
                count = tbl[0]['Count']
                error = tbl[0]['Count_Error']
                All_data.add_row((count,error,final_date))
    
        print >> file , All_data

            #print "For %r \n" %f , All_data
    file.close()
   
    #plot(filters)
    calibration(filters,N)





def plot(filters):
    
    for f in filters:
        
        
        '''Opening our photometry results that we calculated, and star results I just calculated in Some other code'''
        file_plot = ascii.read('photometry_'+str(sn_name)+str(f)+'.txt','r')
        file2_plot = ascii.read('star_'+str(sn_name)+str(f)+'.txt','r')
        
        '''comparing Star results and SN results in one plot for different filter'''
        
        plt.errorbar(file_plot["Time"],file_plot["Count"],file_plot["Error"],file_plot["Error"]*0,fmt='ko',label='SN')
        plt.errorbar(file2_plot["Time"],file2_plot["Count"],file2_plot["Error"],file2_plot["Error"]*0,fmt='o',label='Star')
        plt.legend()
        plt.xlabel("Time-MJD")
        plt.ylabel("Count")
        plt.title(str(sn_name)+str(f))
        plt.savefig(str(sn_name)+str(f))
        plt.show()


def convert_time(analysis):
    
    #analysis is the fits file in the photometry function
    t = Time(analysis[0].header['DATE-OBS'])
    t.format = 'mjd'
    result = t.value
    return result

def calibration(filters,number):
    
    results = open('final_results_'+str(sn_name)+'.csv','w')
    
    for f in filters:
        
        '''Opening our photometry results that we calculated, and star results I just calculated in Some other code'''
        file = ascii.read('photometry_'+str(sn_name)+str(f)+'.txt','r')
        file2 = ascii.read('star_'+str(sn_name)+str(f)+'.txt','r')
        
        
        mag_star=13.018
        
        file["Count_Star"] = file2["Count"]
        file["Inst_Star"] = -2.5*np.log10(file2["Count"])
        file["zero_point"] = mag_star-file["Inst_Star"]
        file["Inst_SN"] = -2.5*np.log10(file["Count"])
        file["Mag_Star"] = file["Inst_Star"]+file["zero_point"]
        file['Mag_sn_'+str(f)] = -2.5*np.log10(file["Count"])+file["zero_point"]
        file["Error_mag"] = 1.086*number/file['Mag_sn_'+str(f)]
        print file
        
        
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


photometry()













