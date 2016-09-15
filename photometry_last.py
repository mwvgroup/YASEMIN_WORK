from astropy.io import fits #For opening fits files
from astropy.coordinates import SkyCoord #For creating sky coordinates
from astropy import units as u #For assigning units to numbers
from photutils import aperture_photometry, SkyCircularAperture #For perfomring photometry
from astropy.cosmology import WMAP9 as cosmo
from astropy.table import Table
import csv
import os.path
import numpy as np

save_path='/Users/zeynepyaseminkalender/Documents/MyFiles_Pyhton/SuperLOTIS_final'
sn_file=np.array(['sn13dadI.fits','sn13dadR.fits','sn13dadV.fits','sn13ddgI.fits','sn13ddgR.fits','sn13ddgV.fits','sn13dgeV.fits','sn13dkjI.fits','sn13dkjR.fits','sn13dkjV.fits','sn13dklI.fits','sn13dklR.fits','sn13dklV.fits','sn13dkxI.fits','sn13dkxR.fits','sn13dkxV.fits','sn13dyR.fits','sn13fjI.fits','sn13fjR.fits','sn13fjV.fits','sn13fnB.fits','sn13fnI.fits','sn13fnR.fits','sn13fnV.fits','snSNhunt206I.fits','snSNhunt206R.fits','snSNhunt206V.fits'])


for x in sn_file:
    
    with fits.open(x) as analysis:
   
        z = .08
        r = 1 * u.kpc / cosmo.kpc_comoving_per_arcmin(z)
    
        cordinate = SkyCoord('01:48:08.66 +37:33:29.22', unit = (u.hourangle, u.deg))
        #We can now create an aperture
        aperture = SkyCircularAperture(cordinate, r)
    
        #Performing photometry returns a table
        phot_table = aperture_photometry(analysis[0], aperture)
        
        # the path is to write the results on upper directory
        completeName=os.path.join(save_path,'photometry.txt')
        file=open(completeName,'w')
        print >> file ,phot_table["aperture_sum"] ,z #taking the residuals and directing them to my file
        file.close()
        print phot_table



