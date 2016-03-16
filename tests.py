import numpy as np
import sys

rf = (1, 3)  #  Reference site`s number " if there is 1 reference site, remove the parenthesis
n = 20  # number of frequencies
fmin = 0.4  # Min frequency
fmax = 15  # Max frequency
bs = 3.5  # S-wave propagation speed in km/h
ps = 2.8  # Density of rocks near source in .......
inputfolder = '/home/babak/PycharmProjects/Simulation' \
                  '/GeneralInversion/Inputs/'  # Path to input folder "/" is needed at the end

outputfolder = '/home/babak/PycharmProjects/Simulation' \
                  '/GeneralInversion/Results/'  # Path to results folder "/" is needed at the end


# ----------------- Checking for input files to be like we want -------------------
count = 0
try:
    ratio = np.genfromtxt(inputfolder + '/ratio.csv', delimiter=',')
    count += 1
except Exception:
    print 'Please check the input files..ratio.csv cannot be found'
    inputfolder = raw_input('Please give me the right path to input files')

try:
    spec = np.genfromtxt(inputfolder + '/spectrums.csv', delimiter=',')
    count += 1
except Exception:
    print 'Please check the input files..spectrums.csv cannot be found'
    inputfolder = raw_input('Please give me the right path to input files')

try:
    earth = np.genfromtxt(inputfolder + '/earth.csv', delimiter=',')
    count += 1
except Exception:
    print 'Please check the input files..earth.csv cannot be found'
    inputfolder = raw_input('Please give me the right path to input files')

# ----------------- Checking in depth of files -----------------
    #  Checkin if 3 first rows of spectrums and ratio file are the same
check = True
if count == 3:
    for i in range(np.shape(ratio)[1]):
            if (ratio[0, i] == spec[0, i]) and (ratio[2, i] == spec[2, i]) and (ratio[1, i] == spec[1, i]):
                pass
            else:
                sys.exit('Spectrums and ratio in column(first 3 rows) number %d are not equal..'
                     'please fix the problem first' % (i))


# ----------------- Cheking if reference sites are correct or nor ---------
