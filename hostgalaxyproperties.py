# This Code is Exercise for Opening Fit files downloaded from Vizier
# and using their data to make plot
# Paper is:
# HOST GALAXY PROPERTIES AND HUBB;E RESIDUALS OF TYPE IA SUPERNOVA FROM THE NEARBY SUPERNOVA FACTORY

from astropy.io import ascii
from astropy.io import fits
from astropy.table import table
import numpy as  np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits

#t=fits.open('hostgalaxyproperties.fit')
#t.info()
t=Table.read('hostgalaxyproperties.fit')
print t
print t["E_B-V_"]

print t.colnames
print t["SN"]
print t["E_B-V_"]
plt.scatter(t["E_B-V_"],t["z"])
plt.xlabel("E_B_V magnitude")
plt.ylabel("red shift")
plt.title("random plot for host galaxies")
plt.show()



