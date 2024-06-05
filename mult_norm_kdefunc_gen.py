## Generates kernel density functions

# For operation with find_ts_loadin_norm_kde.py
# Run this program first
from math import sqrt, pi, exp
from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table,Column
import time
import sys

## Lines to run in terminal
# nohup python3.8 mult_norm_kdefunc_gen.py 080916C 0 1000 0 1000 gaussian gaussian_revision > 080916C_gaussian.out &
# disown ######
## This will allow you to close the terminal while the program is running and it will continue to run while printing the output to 080916C_gaussian.out
## Be sure to create a folder in this directory labeled "kdefuncs" for files to save in

GRB = sys.argv[1]
tmin = float(sys.argv[2])
tmax = float(sys.argv[3])
kdemin = float(sys.argv[4])
kdemax = float(sys.argv[5])
kernel_type = sys.argv[6]
label = sys.argv[7]

# This is the space between points, must be smaller than 0.01 for some reason?
space = 0.0001

data = ascii.read('2-'+GRB+'.txt',format='fixed_width',header_start=5,data_start=6)
#data = ascii.read('random(short)_burst.txt',format='fixed_width',header_start=0,data_start=1)
data1 = data['TIME(Seconds)']
#data1 = data['Time']
data2 = []
data3 = []
kdefunc_values = 0

mult = 1/(kdemax-kdemin)

for l in data1:
    if l<kdemax and l>kdemin:
        withzero = l - kdemin
        adjusted = mult*withzero
        data2.append(adjusted)
        
tmin = (tmin-kdemin)*mult
tmax = (tmax-kdemin)*mult

for u in data2:
    if u<tmax and u>tmin:
        data3.append(u)

#print(data2)
#print(data3)

# Setting the binsizes that we will be creating

#binsizes1 = np.arange(0.0005, 0.1005, 0.0005)
##binsizes1 = np.linspace(1e-5, 5e-5, 41)
##binsizes1 = [1e-5, 2e-5, 5e-5]
binsizes1 = [1e-9, 2e-9, 5e-9, 1e-8, 2e-8, 5e-8, 1e-7, 2e-7, 5e-7, 1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2]
binsizes = [float("{:.12f}".format(num)) for num in binsizes1]
print(binsizes)

# This just causes the program to not start its calculations until 30 seconds has
# elapsed. This allows the user to disown the program before getting spammed with
# outputs.
##time.sleep(30)

# ===================================================
# kde_pdf and kde_cdf are used for compiling kernel |
# density and distribution estimates.               |
# ===================================================


def kde_pdf(data, kernel_func, bandwidth):
    """Generate kernel density estimator over data."""
    kernels = dict()
    n = len(data)
    for d in data:
        kernels[d] = kernel_func(d, bandwidth)
    def evaluate(x):
        """Evaluate `x` using kernels above."""
        pdfs = list()
        for d in data: pdfs.append(kernels[d](x))
        return(sum(pdfs)/n)
    return(evaluate)

# ============================================
# Gaussian Kernel PDF                        |
# ============================================


def gaussian_pdf(x_i, bandwidth):
    """Return Gaussian kernel density estimator."""
    x_bar  = x_i
    def evaluate(x):
        """Evaluate x."""
        pdf = (np.sqrt(2*np.pi*bandwidth**2)**-1) * np.exp(-((x - x_bar)**2)/(2*bandwidth**2))
        return(pdf)
    return(evaluate)

# ============================================
# Triangular Kernel PDF                      |
# ============================================

def triangular_pdf(x_i, bandwidth):
    """Return triangular kernel density estimator."""
    lowerb = (x_i - bandwidth)
    upperb = (x_i + bandwidth)
    def evaluate(x):
        """Evaluate x."""
        if  x <= lowerb: pdf=0
        elif x > upperb: pdf=0
        else: pdf = ((bandwidth-abs(x-x_i))/bandwidth**2)
        return (pdf)
    return (evaluate)

# ============================================
# Uniform Kernel PDF                      |
# ============================================

def uniform_pdf(x_i, bandwidth):
    """Return uniform kernel density estimator."""
    lowerb = (x_i - bandwidth)
    upperb = (x_i + bandwidth)
    def evaluate(x):
        """Evaluate x."""
        if  x<=lowerb: pdf=0
        elif x>upperb: pdf=0
        else: pdf=(1/(2*bandwidth))
        return(pdf)
    return(evaluate)

# ===============================================
# Logic for producing kernel density estimate   |
# vizualizations.                               |
# ===============================================
#import seaborn as sns
#sns.set(color_codes=True)
#plt.rcParams["figure.figsize"] = (15,10)


vals  = data2
xvals = np.arange(min(vals), max(vals), space)

# Beginning of actually creating different kdefuncs using
# different binsizes

for h in binsizes:
    print(h)
    del kdefunc_values
    kdefunc_values=Table(names=('x','y'))

    if kernel_type == "gaussian":
        dist_1 = kde_pdf(data=vals, kernel_func=gaussian_pdf, bandwidth=h)
    if kernel_type == "uniform":
        dist_1 = kde_pdf(data=vals, kernel_func=uniform_pdf, bandwidth=h)
    if kernel_type == "triangular":
        dist_1 = kde_pdf(data=vals, kernel_func=triangular_pdf, bandwidth=h)
    y1 = [dist_1(i) for i in xvals]
    xset = []
    yset = []
    yset_unnorm = []

    for r in range(len(xvals)):
        if xvals[r]<=tmax and xvals[r]>=tmin:
            xset.append(xvals[r])
            yset_unnorm.append(y1[r])

    #normalizing the y-axis so that I can generate random numbers with random.random()
    #instead of random.random()*ymax, which will make my code faster
    for i in yset_unnorm:
        norm_y = i/(max(yset_unnorm))
        yset.append(norm_y)

    for f in range(len(xset)):
        kdefunc_values.add_row([xset[f],yset[f]])
    
    #print(kdefunc_values)
    #print("Number of Photons: "+str(len(data3)))
    print('kdefuncs/'+label+'_'+str(GRB)+'_norm_kdefunc_binsize_'+str(h)+'.txt')

    kdefunc_values.write('kdefuncs/'+label+'_'+str(GRB)+'_norm_kdefunc_binsize_'+str(h)+'.txt',format='ascii.fixed_width',overwrite=True)