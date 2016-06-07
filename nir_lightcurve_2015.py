# NEAR INFRARED LIGHTCURVE FRIEDMAN PAPER 9 TABLE

from astropy.io import ascii
import numpy as  np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits



t1=Table.read('asu1.fit')
t2=Table.read('asu-2.fit')
t3=Table.read('asu-3.fit')
t5=Table.read('asu-5.fit')
t8=Table.read('asu-8.fit')
t9=Table.read('asu-9.fit')
print t1
#print t1["z"]
print t2
print t3
print t5
print t8
print t9

