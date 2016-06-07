# THIS CODE IS FOR USING THE DATA WRITTEN IN TXT FILE
# READING ITS COLUMNS
# DRAWING VARIOUS TYPES OF PLOTS

from astropy.io import ascii
import numpy as numpy
import matplotlib.pyplot as plt
from astropy.table import Table


tbl=ascii.read("host_galaxy.txt")
print ("SURVEY OF SN IA BIAS RESULTS")
print tbl
print tbl["HubbleR"]
print tbl["Median"]
#plt.scatter(tbl["Median"],tbl["HubbleR"])
plt.errorbar(tbl["Median"],tbl["HubbleR"],tbl["Error"],fmt='o')
plt.ylabel("Hubble Residual Magnitude")
plt.xlabel("Median(log(M*/M_star)")
plt.title("Hubble Residul Mag Plot")
plt.show()


