from __future__ import division
import numpy as np
import matplotlib.pyplot as plt



def geometrical(distlist):

    r1 = 0.0
    p1 = -0.2
    r2 = 70.0
    p2 = 0.5
    r3 = 150.0
    rr1 = 1.0
    rr2 = 501.0
    nr = 500

    zr = np.zeros(len(distlist), dtype='float')


    counter = 0
    for i in distlist:
        if r1 <= i < r2:
            zr[counter] = 1.0 / i
        elif r2 <= i < r3:
            zr[counter] = 1.0 / r2*((r2/i)**p1)
        elif i >= r3:
            zr[counter] = 1.0 / r2*((r2/r3)**p1)*((r3/i)**p2)
        counter += 1
    return zr  #, r1, r2, r3


def geometricalplotting(zr):

    fig2 = plt.figure()
    geometricalspearding = plt.loglog(np.arange(0, len(zr)), zr)
    #fig2.savefig('/home/babak/PycharmProjects/Simulation/General-Inversion/Test/'
    #             'geo.eps', format='eps', dpi=1000, facecolor='red', edgecolor='blue')
    plt.show()

# -------------test----------------------
# a = np.genfromtxt('/home/babak/PycharmProjects/Simulation/GeneralInversion/Inputs/spectrums.csv',
#                   dtype='float', delimiter=',')
# print a[2, :]
# zr = geometrical(a[2, :])
# print np.log(zr)
# geometricalplotting(zr)
# -------------test----------------------

