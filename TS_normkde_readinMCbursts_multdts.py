# For operation with mult_kdefunc_gen.py
# Run this program second
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
from math import e, sqrt, pi, exp, log
import math
import time

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
filename = GRB+'_1000_MCbursts_uniformrevision1.txt'

random.seed(seednumber)
rng_state = np.random.get_state()

# TIME DIFFS - CHANGE WHEN PUTTING IN NEW BURST
##minvals = [1e-9, 2e-9, 5e-9, 1e-8, 2e-8, 5e-8, 1e-7, 2e-7, 5e-7, 1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2]
##minvals = [0.01, 0.02, 0.05]
##minvals = [0.05]

# BINSIZES - CHANGE FOR EACH RUN
#binsizes = np.linspace(1e-5, 5e-5, 41)
##binsizes = [1e-5, 2e-5, 5e-5]
##binsizes = [5e-5]
binsizes = [1e-9, 2e-9, 5e-9, 1e-8, 2e-8, 5e-8, 1e-7, 2e-7, 5e-7, 1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2]

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

# This just causes the program to not start its calculations until 30 seconds has
# elapsed. This allows the user to disown the program before getting spammed with
# outputs.
##time.sleep(30)

# Formatting the array of time differences and binsizes
#------------------------------------------------------------------------
#minval_array = [float("{:.12f}".format(num)) for num in minvals]
binsize_array = [float("{:.12f}".format(num)) for num in binsizes]
#------------------------------------------------------------------------

#SETTING UP FINAL DATA TABLE
#------------------------------------------------------------------------
final_results=Table(names=('Binsize','%','Significance'))
#------------------------------------------------------------------------

#LOADING IN MC BURSTS
#------------------------------------------------------------------------
def read_list_from_file(file):
    text = file.read()
    nested_list = eval(text)
    return nested_list

#with open(GRB+'_'+str(trials)+'_MCbursts1.txt', 'r') as file:
with open(filename, 'r') as file:
    MC_burst_list = read_list_from_file(file)
#------------------------------------------------------------------------

# QUANTIFY BUNCHED PHOTONS IN ORIGINAL BURST
#------------------------------------------------------------------------
timediffs = []
for i in range(len(data12)-1):
    dt = data12[i+1] - data12[i]
    timediffs.append(dt)

def take_logarithm_of_array_base_e(array):
    return [math.log(x) / math.log(math.e) for x in array]

logs_of_timediffs = take_logarithm_of_array_base_e(timediffs)

##original_product = math.exp(sum(logs_of_timediffs))
original_product = sum(logs_of_timediffs)
##print('Original Product: '+str(original_product))
#------------------------------------------------------------------------

for h in range(len(binsize_array)):
    ##np.random.set_state(rng_state)
    binsize = binsize_array[h]
    print(binsize)

    #Reset/create m_counter
    m_counter = 0

    #For loop in place of parallelization
    for a in range(trials):

        # QUANTIFY BUNCHED PHOTONS IN MONTECARLO BURST
        #------------------------------------------------------------------------
        xarray = MC_burst_list[h][a]

        timediffs = []
        for i in range(len(xarray)-1):
            dt = xarray[i+1] - xarray[i]
            timediffs.append(dt)

        logs_of_timediffs = take_logarithm_of_array_base_e(timediffs)

        ##MC_product = math.exp(sum(logs_of_timediffs))
        MC_product = sum(logs_of_timediffs)
        ##print('MC Product: '+str(MC_product))
        #------------------------------------------------------------------------

        #COMPARE MC BURST TO ORIGINAL BURST
        #------------------------------------------------------------------------
        if MC_product <= original_product:
            m_counter += 1
        
    binsize_result = m_counter/trials*100
    #------------------------------------------------------------------------
    
    #FIND HOW SIGNIFCANT RESULTS ARE & ADD RESULTS TO TABLE
    #------------------------------------------------------------------------

    answer = 0
    if binsize_result < (((1)/(21.98))*100):
        answer = 2
    if binsize_result < (((1)/(370.))*100):
        answer = 3
    if binsize_result < (((1)/(1744000))*100):
        answer = 5

    print(binsize_result)
    final_results.add_row([binsize,binsize_result,answer])
    print("Time to run: "+str(datetime.now() - startTime))
    #------------------------------------------------------------------------

print(final_results)

#SAVING RESULTS TABLE
#------------------------------------------------------------------------
final_results.write(str(label)+'_multdts_'+str(GRB)+'_timescale_'+str(trials)+'_seed'+str(seednumber)+'.txt',format='ascii.fixed_width',overwrite=True)
#------------------------------------------------------------------------

print("Time to run: "+str(datetime.now() - startTime))