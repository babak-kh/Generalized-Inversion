import numpy as np
import effectscalculations as effects
from interpolating import interpolating as interp
from GeneralizedInversion import gmatrix, svdpart
from  rocksiteamplification import BooreGenericRockSite as bgrs
import os

# #  ------------------- Variables -----------------------------
#
kh = 0.05
kv = 0.025
rf = (77, 16)  # Reference site`s number " if there is 1 reference site, remove the parenthesis
n = 20  # number of frequencies
fmin = 0.4  # Min frequency
fmax = 15  # Max frequency
bs = 3.5  # S-wave propagation speed in km/h
ps = 2.8  # Density of rocks near source in .......
inputfolder = '/home/babak/PycharmProjects/Simulation' \
              '/GeneralInversion/Inputs/'  # Path to input folder "/" is needed at the end

outputfolder = '/home/babak/PycharmProjects/Simulation' \
               '/GeneralInversion/Results/'  # Path to results folder "/" is needed at the end

# --------------------- Making Needed folders ---------------------

if not os.path.exists(outputfolder):  # Making main reults folder
    os.makedirs(outputfolder)
if not os.path.exists(outputfolder + 'Interpolating'):  # Making interpolating folder
    os.makedirs(outputfolder + 'Interpolating')
if not os.path.exists(outputfolder + 'Geometrical-spreading'):  # Making Geometrical spreading folder
    os.makedirs(outputfolder + 'Geometrical-spreading')
if not os.path.exists(outputfolder + 'matrices'):  # Making matrices folder
    os.makedirs(outputfolder + 'matrices')
if not os.path.exists(outputfolder + 'svd'):  # Making SVD folder
    os.makedirs(outputfolder + 'svd')
if not os.path.exists(outputfolder + 'svd/Covariance'):  # Making Covariance folder
    os.makedirs(outputfolder + 'svd/Covariance')
if not os.path.exists(outputfolder + 'Path-effect'):  # Making Path results folder
    os.makedirs(outputfolder + 'Path-effect')
if not os.path.exists(outputfolder + 'Siteamplifications'):  # Making site amplification results folder
    os.makedirs(outputfolder + 'Siteamplifications')
if not os.path.exists(outputfolder + 'Source-effects'):  # maikng source effect results folder
    os.makedirs(outputfolder + 'Source-effects')
if not os.path.exists(outputfolder + 'Source-effects/grid-search'):  # Grid search results folder
    os.makedirs(outputfolder + 'Source-effects/grid-search')

# ------------------- Parts of the program to run ------------------------------

interpolation = False
reference_site_amplification = False
matrices_making = False   # interpolation and reference_site_amplification should be "TRUE"
svd_calculation = False   # interpolation and reference_site_amplification and matrices_making should be "TRUE"
path_effect = False       # interpolation and reference_site_amplification and matrices_making should be "TRUE"
site_effect = False
source_effect = False
grid_search = True

# ------------------ Flow of the program --------------------
# # Calculating reference site amplification using rocksiteamplification code written in Codes foder
if reference_site_amplification:
    referenceamplification = bgrs(fmin, fmax, n)
# # Calling interpolate to make new inputs for spectrum and ratio and frequencies
if interpolation:
    freqs, spec, rat = interp(inputfolder, outputfolder, n, fmin, fmax)
# # Calling gmatrix program in order to make the matrices for svd calculations
if matrices_making:
    g, datamatrix, eqcount, stcount = gmatrix(freqs, spec, rat, outputfolder, rf, referenceamplification, kh, bs, n)
    print 'Please check the number of Earthquakes and Stations'
    print 'Number of earthquakes is %d' % (eqcount)
    print 'Number of stations is %d' % (stcount)
# # Calling svdpart in order to calculate the equation via SVD method
if svd_calculation:
    results = svdpart(outputfolder, g, datamatrix)
    np.savetxt(outputfolder + 'svd/answer.txt', results,
               fmt='%0.5f', delimiter=',')
# # Calculating PATH effects
if path_effect:
    effects.pathcalculations(outputfolder, n)
# #
# # Calculating SITE effects
# # If u wish to have H/V results ploted on the site amplifications, True HtoV variable
if site_effect:
    dk = kh - kv
    HtoV=(False,dk)
    effects.siteextraction(inputfolder, outputfolder, eqcount, stcount, n, HtoV)
# # Calculating SOURCE effects
if source_effect:
    effects.sourceextraction(inputfolder, outputfolder, eqcount, bs, ps, n)
# Grid search method for calculating Sigma, Gamma, Mw
# Range of these variables:
if grid_search:
    sigma = (1, 600, 600)   #  (Start, End, Number of samples)
    gamma = (2.0, 2.0, 1)    #  (Start, End, Number of samples)
    magnitudes = (0.5, 30)   #  (Magnitude increase limit, Step)
    effects.gridsearch(inputfolder, outputfolder, sigma, gamma, magnitudes, bs, n)
##Plotting the results of grid search and source moment rates
# effects.gridsearchplots(inputfolder, outputfolder)
