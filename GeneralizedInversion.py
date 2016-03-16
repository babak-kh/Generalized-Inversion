#!/usr/bin/env python

from numpy.linalg import linalg
import numpy as np
import interpolating
import GeometricalSpreading
# ------------------basic variables---------------------------


# bs = 3.5
# als = 5.98
# muu = 4.41*10**4
# f1 = .4
# f2 = 15
# n = 20
# ps = 2.8
# k0 = .05
# lsf1 = .4
# lsf2 = 2.5
# hsf1 = 3
# hsf2 = 6.3
# c = 2*0.55*1/mt.sqrt(2)/(4*mt.pi*bs**3*ps)*10**(-20)
# a1 = (mt.log(f2)-mt.log(f1))/(n-1)

# -------------- Making G matrix ---------------


def gmatrix(freqs, spec, rat, outputpath, rf, referenceamplification, k0=0.0451, bs=3.5, n=20, f1=0.4, f2=15):

    #print spec
    first3rows = spec[:3, :]
    #print first3rows
    eqcount = np.max(first3rows[0, :])
    stationscount = np.max(first3rows[1, :])

    zr = GeometricalSpreading.geometrical(first3rows[2, :])
    np.savetxt(outputpath + 'Geometrical-spreading/geometrical spreading.txt', np.log(zr),
                fmt='%0.2f', delimiter=',')
    spectrumshape = np.shape(spec)
    recordcount = spectrumshape[1]
    if type(rf) == 'int':
        rfnumber = 1
    else:
        rfnumber = len(rf)
# ----------------------- Making G matrix -------------------------------
    g = np.zeros([recordcount*n+rfnumber*n, (1+eqcount+stationscount)*n],  # len(rf) == 1
                 dtype='float64')  # define the shape of g matrix
    print np.shape(g)
    datamatrix = np.zeros([recordcount*n+len(rf)*n, 1])  # Defines the shape of data matrix
    print np.shape(datamatrix)
    counter = 0
    lengggth = np.shape(first3rows)[1]

# ----------------------------------making g matrix with reference site and path effects --------
    check = True
    for i in range(lengggth):  # i is earthquake number
        eqnumber = first3rows[0, i]
        stnumber = first3rows[1, i]

        for j in range(n):

            sigma = rat[j+3, i]
            # print counter, i, eqnumber, stnumber
            g[counter+j, j] = (spec[2, i] * np.pi * freqs[j] * -sigma) / bs  # Path elements
            g[counter+j, n*(eqnumber-1)+j+n] = sigma  # source elements
            g[counter+j, (eqcount*n)+((stnumber-1)*n)+j+n] = sigma  # site elements

            datamatrix[counter+j, 0] = ((np.log(spec[3+j, i]) + (np.pi * k0 * freqs[j]) - np.log(zr[i]))*sigma)
            if (i == lengggth - 1) and check:
                check = False
                for k in range(rfnumber):  # elements for reference sites in gmatrix and datamatrix
                    if rfnumber == 1:
                        rfk = 1
                    else:
                        rfk = rf[k]
                    for p in range(n):
                        # print k, p
                        datamatrix[(recordcount*n)+p+n*k, 0] = np.log(referenceamplification[p])
                        g[(recordcount*n)+(n*k)+p, (eqcount*n)+((rfk-1)*n)+p+n] = 1
        counter += n
    string = 'Number of reference sites are %d\n' \
             'Data matrix has %d rows and %d columns\n' \
             'G matrix has %d rows and %d columns' % (rfnumber, np.shape(datamatrix)[0], np.shape(datamatrix)[1],
                                                      np.shape(g)[0], np.shape(g)[1])
    np.savetxt(outputpath + 'matrices/gmatrix.txt', g,  fmt='%.4e',
               delimiter=',')
    np.savetxt(outputpath + 'matrices/datamatrix.txt', datamatrix,  fmt='%.4e',
               delimiter=',')
    with open(outputpath + 'matrices/G-Test.txt', "w") as text_file:
        text_file.write(string)
    return g, datamatrix, eqcount, stationscount

# ------------------------- SVD decomposition and source-site-path results -------------------


def svdpart(outputpath, finalgmatrix, finaldatamatrix):

    u, s, v = linalg.svd(finalgmatrix, full_matrices=False)
    np.savetxt(outputpath + 'svd/u.txt', u,  fmt='%.4e',
               delimiter=',')

    np.savetxt(outputpath + 'svd/v.txt', np.transpose(v),  fmt='%.4e',
               delimiter=',')
    for i in range(len(s)):
        s[i] = 1.0 / s[i]
        # print s[i]
        if s[i] > 10 ** 6:
            s[i] = 0
    S = np.diag(s)

    np.savetxt(outputpath + 'svd/s.txt', S,  fmt='%.4e',
               delimiter=',')

    result = np.dot(np.dot(np.dot(np.transpose(v),S), np.transpose(u)), finaldatamatrix)
    # --------------- Covariance matrix calculation ----------------------------
    a = np.dot(np.dot(np.dot(np.transpose(v), S), S), v)
    a2 = np.dot(finalgmatrix, result)
    a3 = (finaldatamatrix - a2)
    ng = np.shape(finalgmatrix)
    sig2d = np.dot(np.transpose(a3), a3/(ng[0]-ng[1]))
    covd = a * sig2d
    np.savetxt(outputpath + 'svd/Covariance/a.txt', a,  fmt='%.2e',
               delimiter=',')
    np.savetxt(outputpath + 'svd/Covariance/a2.txt', a2,  fmt='%.2e',
               delimiter=',')
    np.savetxt(outputpath + 'svd/Covariance/a3.txt', a3,  fmt='%.2e',
               delimiter=',')
    np.savetxt(outputpath + 'svd/Covariance/covd.txt', covd,  fmt='%.2e',
               delimiter=',')
    print sig2d
    return result








