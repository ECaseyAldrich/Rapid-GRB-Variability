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
from math import sqrt, pi, exp
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
minvals = [1e-9, 2e-9, 5e-9, 1e-8, 2e-8, 5e-8, 1e-7, 2e-7, 5e-7, 1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2]
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
minval_array = [float("{:.12f}".format(num)) for num in minvals]
binsize_array = [float("{:.12f}".format(num)) for num in binsizes]
#------------------------------------------------------------------------

arrays = 0 #Giving my code something to delete in 5 lines :)

#SETTING UP FINAL DATA TABLE
#------------------------------------------------------------------------
final_results=Table(names=('Binsize','Lowest %','Significance'))
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


for h in range(len(binsize_array)):
    ##np.random.set_state(rng_state)
    binsize = binsize_array[h]
    print(binsize)

    del arrays

    #Reset/create m_counter
    m_counter = []
    for k in range(len(minval_array)):
        m_counter.append(0)

    #For loop in place of parallelization
    for a in range(trials):

        #Setting up temporary m_counter
        temp_m_counter = []
        for k in range(len(minval_array)):
            ##temp_m_counter = np.append(temp_m_counter, 0)
            temp_m_counter.append(0)

        arrays = []

        for j in range(len(minval_array)):

            minval = minval_array[j]

            # COUNT BUNCHED PHOTONS IN ORIGINAL BURST
            #------------------------------------------------------------------------
            nearx1 = []
            for k in range(len(data12)):
                if k == 0:
                    dt2 = abs(data12[k+1]-data12[k])
                elif k == (len(data12)-1):
                    dt2 = abs(data12[k]-data12[k-1])
                else:
                    dt2 = min(abs(data12[k]-data12[k-1]), abs(data12[k]-data12[k+1]))
                if dt2<minval:
                    nearx1.append(data12[k])

            fraction=((len(nearx1))/(len(data12)))
            #------------------------------------------------------------------------

            #COUNT BUNCHED PHOTONS IN RANDOM BURST
            #------------------------------------------------------------------------
            nearx, neary = [], []
            xarray = MC_burst_list[h][a]

            for k in range(len(xarray)):
                if k == 0:
                    dt2 = abs(xarray[k+1]-xarray[k])
                elif k == (len(xarray)-1):
                    dt2 = abs(xarray[k]-xarray[k-1])
                else:
                    dt2 = min(abs(xarray[k]-xarray[k-1]), abs(xarray[k]-xarray[k+1]))
                if dt2<minval:
                    nearx.append(xarray[k])

            fraction2 = ((len(nearx))/len(xarray))
            #------------------------------------------------------------------------

            #COMPARE MC BURST TO ORIGINAL BURST
            #------------------------------------------------------------------------
            if fraction2 >= fraction:
                temp_m_counter[j] += 1

        ##m_counter = m_counter + temp_m_counter
        #print(m_counter, temp_m_counter)
        ##print("      "+str(m_counter))

        for k in range(len(minval_array)):
            m_counter[k] += temp_m_counter[k]
            
        ##print("      "+str(m_counter))
            
        arrays.append(m_counter)
        
    binsize_results = [x / trials*100 for x in arrays[-1]]
    #------------------------------------------------------------------------
    
    #FIND HOW SIGNIFCANT RESULTS ARE & ADD RESULTS TO TABLE
    #------------------------------------------------------------------------
    ##k = min(binsize_results)
    ##Finding the second smallest value in the array.
    def find_second_smallest(arr):
        sorted_arr = sorted(arr)
        return sorted_arr[1]
    
    k = find_second_smallest(binsize_results)

    answer = 0
    if k < (((1)/(21.98))*100):
        answer = 2
    if k < (((1)/(370.))*100):
        answer = 3
    if k < (((1)/(1744000))*100):
        answer = 5

    print(binsize_results)
    final_results.add_row([binsize,k,answer])
    print("Time to run: "+str(datetime.now() - startTime))
    #------------------------------------------------------------------------

print(final_results)

#SAVING RESULTS TABLE
#------------------------------------------------------------------------
final_results.write(str(label)+'_'+str(GRB)+'_timescale_'+str(trials)+'_seed'+str(seednumber)+'.txt',format='ascii.fixed_width',overwrite=True)
#------------------------------------------------------------------------

print("Time to run: "+str(datetime.now() - startTime))