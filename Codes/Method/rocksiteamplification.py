import sys
import numpy as np


def BooreGenericRockSite(f1=0.4, f2=15, n=20):


    a1 = (np.log(f2)-np.log(f1))/(n-1)
    w1 = np.arange(0, n, dtype=float)
    w1[0] = f1
    for i in range(1, int(n)):
        w1[i] = np.exp(np.log(w1[0])+(i*a1))

    frequency = np.array([0.01,0.09,0.16,0.51,0.84,1.25,2.26,3.17,6.05,16.6,61.2], dtype=float)
    amplification = np.array([1.00,1.10,1.18,1.42,1.58,1.74,2.06,2.25,2.58,3.13,4.00], dtype=float)
    z = np.polyfit(np.log(frequency), np.log(amplification), 1)
    a = np.poly1d(z)
    rocksiteamplification = [np.exp(a(x)) for x in np.log(w1)]

    log = 'Reference site amplifications are following Boore`s generic rock site amplification.' \
          'Check that out for more information'


    return rocksiteamplification, log

#BooreGenericRockSite()