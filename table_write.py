# This Code is creating a data with random number generator.
#Writing Tables, Reading Tables from dat extention files...

from astropy.io import ascii
import numpy as np
from astropy.table import Table
x=np.random.random(10)
y=x**2
data=Table()
ascii.write([x,y],'values.dat',names=['x','y'])

file=open('values.dat','r')
print file.read()


