def A(kh, kv, fmin, fmax, bs, ps, rf, n,alphas, inputpath, outputpath, interpolationstatus,
    pgastatus, rfamplificationstatus, matricesstatus, svdstatus, pathstatus,
    sitestatus, sourcestatus, gristatus, gridplotstatus, extraplotsstatus):
    import numpy as np
    import effectscalculations as effects
    from interpolating import interpolating as interp
    from GeneralizedInversion import gmatrix, svdpart
    from  rocksiteamplification import BooreGenericRockSite as bgrs
    import os

    import GeneralInversion.Codes.pyqttraining

    # #  ------------------- Variables -----------------------------
    #
    kh = float(kh)  # Kappa for horizontal components
    kv = float(kv)  # Kappa for vertical components
    rf = (16, 77)  # Reference site`s number " if there is 1 reference site, remove the parenthesis
    n = int(20)  # number of frequencies
    fmin = float(fmin)  # Min frequency
    fmax = float(fmax)  # Max frequency
    bs = float(bs)  # S-wave propagation speed in km/h
    alphas = float(alphas) #P-wave propagation velosity in km/h
    ps = float(ps)  # Density of rocks near source in .......
    inputfolder = str(inputpath)  # Path to input folder "/" is needed at the end

    outputfolder = str(outputpath)  # Path to results folder "/" is needed at the end

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
    if not os.path.exists(outputfolder + 'Source-effects/Extra-plots'):  # Grid search results folder
        os.makedirs(outputfolder + 'Source-effects/Extra-plots')
    # ------------------- Parts of the program to run ------------------------------

    interpolation = True
    reference_site_amplification = True
    pga = True
    matrices_making = True   # interpolation and reference_site_amplification should be "TRUE"
    svd_calculation = False   # interpolation and reference_site_amplification and matrices_making should be "TRUE"
    noteffectcalculations = True  # Whole effects calculation process can be ignored with True ing this variable
    if noteffectcalculations:
        path_effect = False     # interpolation and reference_site_amplification and matrices_making should be "TRUE"
        site_effect = False     # Same as Above
        source_effect = False
        grid_search = False
        grid_search_plotting = False
        extrasourceplottings = True
    else:
        path_effect = True     # interpolation and reference_site_amplification and matrices_making should be "TRUE"
        site_effect = True     # Same as Above
        source_effect = True
        grid_search = True
        grid_search_plotting = True
        extrasourceplottings = True
    # ------------------ Flow of the program --------------------
    # # Calculating reference site amplification using rocksiteamplification code written in Codes foder
    if rfamplificationstatus:
        referenceamplification = bgrs(fmin, fmax, n)

    # # Calling interpolate to make new inputs for spectrum and ratio and frequencies
    if interpolationstatus:
        freqs, spec, rat = interp(inputfolder, outputfolder, n, fmin, fmax)
    # # Calculating PGA X Epicentral distance plots
    if pgastatus:
        effects.PGA(outputfolder, spec)
    # # Calling gmatrix program in order to make the matrices for svd calculations
    if matricesstatus:
        g, datamatrix, eqcount, stcount = gmatrix(freqs, spec, rat, outputfolder, rf, referenceamplification, kh, bs, n)
        print 'Please check the number of Earthquakes and Stations'
        print 'Number of earthquakes is %d' % (eqcount)
        print 'Number of stations is %d' % (stcount)
    # # Calling svdpart in order to calculate the equation via SVD method
    if svdstatus:
        results = svdpart(outputfolder, g, datamatrix)
        np.savetxt(outputfolder + 'svd/answer.txt', results,
                   fmt='%0.5f', delimiter=',')
    # # Calculating PATH effects
    if pathstatus:
        effects.pathcalculations(outputfolder, n)
    # #
    # # Calculating SITE effects
    # # If u wish to have H/V results ploted on the site amplifications, True HtoV variable
    if sitestatus:
        dk = kh - kv
        HtoV = (False, dk)
        effects.siteextraction(inputfolder, outputfolder, eqcount, stcount, n, HtoV)
    # # Calculating SOURCE effects
    if sourcestatus:
        effects.sourceextraction(inputfolder, outputfolder, eqcount, bs, ps, n)
    # Grid search method for calculating Sigma, Gamma, Mw
    # Range of these variables:
    if gristatus:
        sigma = (1, 600, 600)   #  (Start, End, Number of samples)
        gamma = (2.0, 2.0, 1)    #  (Start, End, Number of samples)
        magnitudes = (0.5, 40)   #  (Magnitude increase limit, Step)
        effects.gridsearch(inputfolder, outputfolder, sigma, gamma, magnitudes, bs, n, ps, alphas)
    ##Plotting the results of grid search and source moment rates
    if gridplotstatus:
        effects.gridsearchplots(inputfolder, outputfolder, eqcount, stcount)

    ## Extra source plots -------------
    if extraplotsstatus:
        effects.extrasourceplots(inputfolder, outputfolder, eqcount, stcount, n, ps, bs, alphas)

if __name__ == '__main__':
    A()