## Take the kdefuncs, and generates the MC bursts

import sys
from astropy.io import ascii
from astropy.table import Table,Column
from astropy.io import fits
import numpy as np
import astropy.units as u
import random
from random import uniform
#from math import *
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import interp1d
from scipy import stats
import statsmodels.nonparametric.api as smnp
from datetime import datetime
from math import sqrt, pi, exp
import time
from tqdm import tqdm

## Change any variables within this code. Be sure to use the same "label" that you used while generating the corresponding kdefuncs
## This will save all MC bursts to 1 text file that can be read in to "TS_normkde_readinMCbursts.py" and "TS_normkde_readinMCbursts_multdts.py"
## Run these lines in terminal:
##	nohup python3.8 norm_random_number_gen.py > 080916C_burstgeneration.out &
##	disown ######

startTime = datetime.now()

#By Photon Method

GRB = '160709A'

##random.seed(1)

tmin = 0
tmax = 1000
kdemin = 0
kdemax = 1000
trials = 1000
label = 'uniformrevision'
mult = 1/(kdemax-kdemin)
seednumber = 1

random.seed(seednumber)
rng_state = np.random.get_state()


#Reading in GRB photon list & trims the list to the desired time range
data2 = ascii.read('2-'+GRB+'.txt',format='fixed_width',header_start=5,data_start=6)    #define GRB data as table
data1=[]
data12=[]

for row in data2:                                                                       #filter out photons outside of the desired time range
    photontime=row['TIME(Seconds)']
    #photontime = row['Time']
    if photontime<=kdemax and photontime>=kdemin:
        data1.append(photontime)
        
for y in data1:
    if y<=tmax and y>=tmin:
        data12.append(y*mult)

# TIME DIFFS - CHANGE WHEN PUTTING IN NEW BURST
#minvals = [1e-9, 2e-9, 5e-9, 1e-8, 2e-8, 5e-8, 1e-7, 2e-7, 5e-7, 1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2]
##minvals = [0.01, 0.02, 0.05]
##minvals = [0.05]

# BINSIZES - CHANGE FOR EACH RUN
#binsizes = np.linspace(1e-5, 5e-5, 41)
##binsizes = [1e-5, 2e-5, 5e-5]
##binsizes = [5e-5]
binsizes = [1e-9, 2e-9, 5e-9, 1e-8, 2e-8, 5e-8, 1e-7, 2e-7, 5e-7, 1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2]

#binsizes = [1e-9, 2e-9, 5e-9, 1e-8, 2e-8, 5e-8, 1e-7, 2e-7]
#binsizes = [5e-7, 1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4]
#binsizes = [2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2]

# Formatting the array of time differences and binsizes
#------------------------------------------------------------------------
#minval_array = [float("{:.12f}".format(num)) for num in minvals]
binsize_array = [float("{:.12f}".format(num)) for num in binsizes]
#------------------------------------------------------------------------

#GENERATE MC BURSTS
#------------------------------------------------------------------------

MCburstlist = []

for i in range(len(binsizes)):
    print(binsizes[i])
    MCburstlist.append([])
    h = binsize_array[i]

    #READ IN KERNEL DENSITY FUNCTION
    #------------------------------------------------------------------------
    kdedata = ascii.read('kdefuncs/'+label+'_'+str(GRB)+'_norm_kdefunc_binsize_'+str(h)+'.txt',format='fixed_width',header_start=0,data_start=1)
    #print(kdedata)
    kdex = kdedata['x']
    kdey = kdedata['y']
    #------------------------------------------------------------------------

# GENERATE ALL MC BURSTS FOR THIS T_SMOOTH
    for j in tqdm(range(trials)):
        counter = 0
        xarray = []
        yarray = []

        valid_coord = []
        f = interp1d(kdex, kdey, bounds_error=False, fill_value="extrapolate")
        while len(valid_coord) < len(data12):
            x = np.random.uniform(0, 1)
            y = np.random.uniform(0, 1)
            #*ymax
            counter +=1 
            if y < f(x):
                valid_coord.append((x, y))

        valid_coord.sort()

        for k in range(len(valid_coord)):
            xarray.append(valid_coord[k][0])
            yarray.append(valid_coord[k][1])

        ##print(xarray)
        ##print(j)
        MCburstlist[i].append(xarray)
        #print(counter)

#print(MCburstlist[0][0])

def write_list_to_file(my_list, file):
    for item in my_list:
        if isinstance(item, list):
            file.write("[")
            write_list_to_file(item, file)
            file.write("]")
        else:
            file.write(str(item))
        file.write(", ")

with open(GRB+'_'+str(trials)+'_MCbursts_'+label+'1.txt', 'w') as file:
    write_list_to_file(MCburstlist, file)

print("Time to run: "+str(datetime.now() - startTime))