from astropy.io import fits #For opening fits files
from astropy.coordinates import SkyCoord #For creating sky coordinates
from astropy import units as u #For assigning units to numbers
from photutils import aperture_photometry, SkyCircularAperture , CircularAperture#For perfomring photometry
from astropy.cosmology import WMAP9 as cosmo
from astropy.table import Table
from astropy.io import ascii
from glob import glob
import os.path
import os
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from pylab import *
from matplotlib.ticker import ScalarFormatter,FormatStrFormatter


sn_name=raw_input('enter the supernova name: >')

def extract_date_from_fullpath(name):
    """Extract date from the director name of the supernova name.
        
        >>> extract_date_from_fullpath('a/b/c')
        'b'
        >>> extract_date_from_fullpath('/Users/zeynepyaseminkalender/Documents/MyFiles_Pyhton/SuperLOTIS_final/130422/sn13dadV.fits')
        '130422'
        """
    date_extract=name.split(os.sep)
    date=date_extract[6]
    return date

def photometry():
    filters=['I','B','V','R'] #filter names
    for f in filters:
        All_data=Table(names=('Count','Error','Time'))
        file=open('photometry_'+str(sn_name)+str(f)+'.txt','w')
        super_lotis_path='/Users/zeynepyaseminkalender/Documents/MyFiles_Pyhton/SuperLOTIS_final'
        date_search_path=os.path.join(super_lotis_path, '13*')
        search_str=os.path.join(date_search_path,str(sn_name)+str(f)+'.fits')
        for name in glob(search_str):
            date=extract_date_from_fullpath(name)
            final_date=convert_time(date)
            with fits.open(name) as analysis:
                z = .086
                r = 1 * u.kpc / cosmo.kpc_comoving_per_arcmin(z)
                cordinate = SkyCoord('01:48:08.66 +37:33:29.22', unit = (u.hourangle, u.deg))
                aperture = SkyCircularAperture(cordinate, r)
                exp_time= analysis[0].header['EXPTIME'] # error calculation
                
                data_error = np.sqrt(analysis[0].data*exp_time) / exp_time
                tbl=Table(aperture_photometry(analysis[0],aperture,error=data_error),names=('Count','Count_Error','x_center','y_center','center_input'))
                tbl.keep_columns(['Count','Count_Error'])
                count=tbl[0]['Count']
                error=tbl[0]['Count_Error']
                All_data.add_row((count,error,final_date))
    
        print >> file , All_data
        file.close()
    plot(filters)

  #calling the functions


def plot(filters):
    for f in filters:
         fig=plt.figure()
         file_plot= ascii.read('photometry_'+str(sn_name)+str(f)+'.txt','r')
         file2_plot=ascii.read('star_'+str(sn_name)+str(f)+'.txt','r')
         plt.figure(str(f))
         plt.errorbar(file_plot["Time"],file_plot["Count"],file_plot["Error"],file_plot["Error"]*0,fmt='ko')
         plt.errorbar(file2_plot["Time"],file2_plot["Count"],file2_plot["Error"],file2_plot["Error"]*0,fmt='o')
         plt.title(str(sn_name)+str(f))
         ax=fig.add_subplot()
         ax.xaxis.set_major_formatter()
         plt.savefig(str(sn_name)+str(f))
         plt.show()



def convert_time(mytime2):
    
    """ This Function is convert the date we have into MJD units ..
        exp: 131009= 2013 October 9
        >>> 2012-10-09 (this is the format that Time function in astopry accepts.
        >>> 50987.00 is the final time
        
    """

    word2='20'
    for i in range(0,len(mytime2),2):
        word2+=str(mytime2[i]+mytime2[i+1])+'-'
    word2=word2[:-1]
    t=Time(word2)
    t.format='mjd'
    result=t.value
    return result

photometry()













