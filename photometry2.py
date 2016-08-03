import numpy as np
import os
from astropy.io import fits
from astropy.cosmology import WMAP9 as cosmo

from astropy.table import Table
import astropy.io.ascii as asc

from astropy import units as u
from astropy.io.fits import getheader
from photutils import aperture_photometry,SkyCircularAperture,CircularAperture
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS




#SkyCircular Apperture

def photometry(fits_file,radius):
        if os.path.isfile(fits_file):
        
            hdulist=fits.open(fits_file)
            ra = hdulist[0].header["RA"]
            dec=hdulist[0].header["DEC"]
            wcs=WCS(fits_file)
            coords = [ra,dec]
            w = wcs.all_world2pix(ra,dec,1)
            r=radius*u.kpc
            aperture=SkyCircularAperture(w,r*u.arcsec)
            exp_time=hdulist[0].header["EXPTIME"]
            photo_error=np.sqrt(hdulist[0].data/exp_time)
            photo_error=np.sqrt(hdulist[0].data/exp_time)
            phot_table=aperture_photometry(hdulist[0],aperture,error=photo_error)
        else:
            print ("it is not a fit file")

    
def main():
    
    photometry('snPSNJ21134481B.fits',9.44)

    
main()