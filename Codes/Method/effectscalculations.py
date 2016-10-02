import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 1.75
mpl.rcParams['font.size'] = 12
mpl.rcParams['font.family'] = 'serif'
import numpy as np
from GeneralInversion.Codes.Method import GeneralizedInversion
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as nppoly
import math as mt
import xlrd

def PGA(outputpath, spec):

    peakaccel = np.amax(spec[3:,:], axis=0)


    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.xaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.1, zorder=0)
    ax.yaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.1, zorder=1)
    ax.scatter(spec[2,:], peakaccel, s=100,
               color='deepskyblue', edgecolor='black', alpha=0.5, zorder=3)
    ax.set_ylim([-1, 50])
    fig.suptitle('Peak ground acceleration')
    ax.xaxis.set_label_text('Distance(Km)', size=18, family='freeserif', color='#000099')
    ax.yaxis.set_label_text('PGA(mg)', size=18, family='freeserif',
                             color='#000099')
    fig.savefig(outputpath + 'Source-effects/Extra-plots/'
                   + 'PGA' + '.pdf', format='pdf', dpi=200)

def gettingresults():  # This function is not used in main flow of the program
    gmatrice, datamatrice = GeneralizedInversion.gmatrix(51, 17, rf=(1, 3))
    np.savetxt('/home/babak/PycharmProjects/Simulation/GeneralInversion/Results/G.txt', gmatrice,
               fmt='%0.5f', delimiter=',')
    np.savetxt('/home/babak/PycharmProjects/Simulation/GeneralInversion/Results/bsv.txt', datamatrice,
               fmt='%0.2f', delimiter=',')
    result = GeneralizedInversion.svdpart(gmatrice, datamatrice)
    np.savetxt('/home/babak/PycharmProjects/Simulation/GeneralInversion/Results/answer.txt', result,
               fmt='%0.5f', delimiter=',')
    return result


def pathcalculations(outputpath, n=20):
    result = np.genfromtxt(outputpath + 'svd/answer.txt',
                           delimiter=',')
    freqs = np.genfromtxt(outputpath + 'Interpolating/freqs.txt',
                          delimiter=',')
    qfactor = result[:n]
    qfactor = 1.0 / qfactor
    np.savetxt(outputpath + 'Path-effect/q.txt', qfactor,
               fmt='%0.5f', delimiter=',')
    fit = nppoly.polyfit(np.log(freqs), np.log(qfactor), 1)
    aa = np.exp(fit[0])

    # ----------------------------- plotting -------------------

    fig = plt.figure(1)
    fig.suptitle('Q-factor', fontsize=14, fontweight='bold', family='freeserif')
    ax1 = fig.add_subplot(111)
    ax1.xaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.1, zorder=0)
    ax1.yaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.1, zorder=1)

    ax1.loglog(freqs, qfactor, color='#202020', linewidth=2,
               linestyle='--', label='Calculated Qs', zorder=4)
    ax1.loglog(freqs, np.exp(nppoly.polyval(np.log(freqs), fit)),
               alpha=0.5, color='#994c00', linewidth=2,
               label=r'Fitted line  $Q_{s}=%0.2f \times f^{%0.1f}$' % (aa, fit[1]), zorder=3)
    # --------------- Plot Properties ------------------
    ax1.xaxis.set_label_text('Frequency(Hz)', size=12, family='freeserif', color='#000099')
    ax1.yaxis.set_label_text('Q-factor(Qs)', size=12, family='freeserif', color='#000099')

    ax1.get_xaxis().tick_bottom()
    ax1.get_yaxis().tick_left()

    # ax1.spines['top'].set_visible(False)
    # ax1.spines['right'].set_visible(False)

    ax1.set_xlim([0.2, 70])
    ax1.set_ylim([20, 2000])
    ax1.legend(loc='upper left', frameon=False, fancybox=True, shadow=True,
               prop={'size': 14, 'family': 'freeserif'})
    # ------------------- End --------------------------
    fig.savefig(outputpath + 'Path-effect/Quality-factor.pdf', format='pdf', dpi=1000)
    plt.close(fig)
    print '%0.2ff^%0.2f' % (np.exp(fit[0]), fit[1])

    # -------------------------- end of plotting ----------


def siteextraction(inputpath, outputpath, eqcount, stcount, n=20, HtoV=(False, 0)):
    #  Calculating H/V if the variable is True

    w2 = xlrd.open_workbook(inputpath +
                            '/w2.xls')
    sheet = w2.sheet_by_index(0)
    Freq = sheet.col_values(0)
    Freq = np.asarray(Freq)
    stations = np.genfromtxt(inputpath + 'stations final.csv',
                             delimiter=',', dtype=str)
    result = np.genfromtxt(outputpath + 'svd/answer.txt',
                           delimiter=',')
    freqs = np.genfromtxt(outputpath + 'Interpolating/freqs.txt',
                          delimiter=',')
    covd = np.genfromtxt(outputpath + 'svd/Covariance/covd.txt', defaultfmt='.4e', delimiter=',')
    stcount = int(stcount)
    eqcount = int(eqcount)
    # print HtoV
    if HtoV[0]:

        H = np.genfromtxt(inputpath + 'spectrums.csv',
                          delimiter=',')
        V = np.genfromtxt(inputpath + 'Vspectrums.csv',
                          delimiter=',')
        htov = (H[2:, :] * 1.0) / ((V[2:, :] / (2 * np.pi)) * 1.0)
        np.savetxt(outputpath + 'Siteamplifications/H-v.txt',
                   htov, fmt='%0.2f', delimiter=',')
        for i in range(1, stcount + 1):
            counterhtov = 0
            sum = np.zeros(np.shape(H)[0] - 2)
            for j in range(np.shape(H)[1]):
                if V[1, j] == i:
                    sum[:] += htov[:, j]
                    counterhtov += 1
            sum[0] = counterhtov
            sum[:] /= sum[0]
            for ww in range(n):
                sum[ww + 1] *= np.exp(np.pi * (HtoV[1]) * Freq[ww])
            # print sum
            H[2:, i - 1] = sum[:]
        HtoVinterp = np.zeros((n + 3, np.shape(H)[1]))
        for i in range(0, np.shape(H)[1]):
            HtoVinterp[3:, i] = np.interp(freqs, Freq, H[3:, i])
        np.savetxt(outputpath + 'Siteamplifications/HtoVinterp.txt',
                   HtoVinterp, fmt='%0.5f', delimiter=',')
        np.savetxt(outputpath + 'Siteamplifications/HtoVnotinterp.txt',
                   H, fmt='%0.5f', delimiter=',')

    siteresults = np.zeros([n, stcount])
    siteresultspositive = np.zeros([n, stcount])
    siteresultsnegative = np.zeros([n, stcount])
    # print n
    for i in range(stcount):
        for j in range(n):
            siteresults[j, i] = np.exp(result[n + (eqcount * n) + (i * n) + j])
            siteresultspositive[j, i] = np.exp(result[n + (eqcount * n) + (i * n) + j] + \
                                               mt.sqrt(covd[n + (eqcount * n) + (i * n) + j, n + (eqcount * n) + (
                                               i * n) + j]))
            siteresultsnegative[j, i] = np.exp(result[n + (eqcount * n) + (i * n) + j] - \
                                               mt.sqrt(covd[n + (eqcount * n) + (i * n) + j, n + (eqcount * n) + (
                                               i * n) + j]))
    np.savetxt(outputpath + 'Siteamplifications/siteresults.txt',
               siteresults, fmt='%0.5f', delimiter=',')
    np.savetxt(outputpath + 'Siteamplifications/siteresultspstd.txt',
               siteresultspositive, fmt='%0.5f', delimiter=',')
    np.savetxt(outputpath + 'Siteamplifications/siteresultsnstd.txt',
               siteresultsnegative, fmt='%0.5f', delimiter=',')
    counter = 0
    for i in range(abs(stcount) / 9 + 1):
        fig = plt.figure(i)
        # fig.suptitle('Site Amplifications', fontsize=14, fontweight='bold', family='freeserif')
        for j in range(9):
            if counter >= stcount:
                break
            ax1 = fig.add_subplot(3, 3, j + 1)
            if j == 0 or j == 3 or j == 6:
                # ax1.xaxis.set_label_text('Frequency(Hz)', size=12,
                #                          family='freeserif', color='#000099')
                ax1.yaxis.set_label_text('Amplification', size=12,
                                         family='freeserif', color='#000099', )
            ax1.xaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.1, zorder=0)
            ax1.yaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.1, zorder=1)
            if HtoV[0]:
                ax1.loglog(freqs, HtoVinterp[3:, counter], lw=1, ls='-', color='Blue')
            ax1.loglog(freqs, siteresults[:, counter], lw=2, ls='-', color='#FF8000')
            ax1.loglog(freqs, siteresultsnegative[:, counter], lw=2, ls='--',
                       zorder=2, color='#193300')
            ax1.loglog(freqs, siteresultspositive[:, counter], lw=2, ls='--',
                       zorder=3, color='#193300')
            ax1.set_title(str(stations[counter, 0]))
            ax1.set_xlim([0.1, 100])
            ax1.set_ylim([0.1, 30])
            ax1.yaxis.set_ticks([10, 20])
            ax1.yaxis.set_ticklabels([10, 20])
            counter += 1
            # print i
            # print counter
            # print str(i)
        plt.tight_layout()
        fig.savefig(outputpath + 'Siteamplifications/'
                    + str(i) + '.pdf', dpi=200)
        plt.close(fig)

        # plt.show()


def sourceextraction(inputpath, outputpath, eqcount, bs=3.5, ps=2.8, n=20):
    c = 2 * 0.55 * 1 / np.sqrt(2) / (4 * np.pi * bs ** 3 * ps) * 10 ** (-20)
    eqcount = int(eqcount)
    earthquakes = np.genfromtxt(inputpath + 'earth.csv',
                                delimiter=',')

    result = np.genfromtxt(outputpath + 'svd/answer.txt',
                           delimiter=',')
    freqs = np.genfromtxt(outputpath + 'Interpolating/freqs.txt',
                          delimiter=',')
    covd = np.genfromtxt(outputpath + 'svd/Covariance/covd.txt',
                         defaultfmt='.4e', delimiter=',')

    magnitudes = np.zeros(eqcount)
    momentrate = np.zeros([n, eqcount])
    momentratenstd = np.zeros([n, eqcount])
    momentratepstd = np.zeros([n, eqcount])
    efromsvd = np.zeros([n, eqcount])

    for i in range(eqcount):
        for j in range(n):
            efromsvd[j, i] = np.exp(result[n + (i * n) + j])
            momentrate[j, i] = np.exp(result[n + (i * n) + j]) / (4 * c * (np.pi ** 2) * (freqs[j] ** 2))
            #print np.exp(result[n + (50 * n) + 19]) / (4 * c * (np.pi ** 2) * (freqs[19] ** 2))
            # print i, j, np.exp(efromsvd[j, i])
            momentratepstd[j, i] = np.exp(
                (result[n + (i * n) + j]) + mt.sqrt(covd[(n + (i * n) + j), (n + (i * n) + j)])) / (
                                   4 * c * (np.pi ** 2) * (freqs[j] ** 2))
            momentratenstd[j, i] = np.exp(
                (result[n + (i * n) + j]) - mt.sqrt(covd[(n + (i * n) + j), (n + (i * n) + j)])) / (
                                   4 * c * (np.pi ** 2) * (freqs[j] ** 2))
            # print np.exp((result[n+(i*n)+j]) + mt.sqrt(covd[(n+(i*n)+j), (n+(i*n)+j)]))
        magnitudes[i] = earthquakes[i, 12]

    moment = magnitudes
    # print moment
    np.savetxt(outputpath + 'Source-effects/efromsvd.txt',
               efromsvd, fmt='%0.5f', delimiter=',')
    np.savetxt(outputpath + 'Source-effects/momentrates.txt',
               momentrate, fmt='%0.5e', delimiter=',')
    np.savetxt(outputpath + 'Source-effects/magnitudesextracted.txt',
               moment, fmt='%0.5f', delimiter=',')
    np.savetxt(outputpath + 'Source-effects/momentratespstd.txt',
               momentratepstd, fmt='%0.5f', delimiter=',')
    np.savetxt(outputpath + 'Source-effects/momentratesnstd.txt',
               momentratenstd, fmt='%0.5f', delimiter=',')
    counter = 0
    for i in range(eqcount / 4 + 1):
        if counter == eqcount:
            break
        fig = plt.figure(i)
        for j in range(4):
            if counter == eqcount:
                break
            ax1 = fig.add_subplot(2, 2, j + 1)
            if j == 0 or j == 2:
                ax1.yaxis.set_label_text(r'Moment rate spectrum $(dyne \times Cm)$', size=12,
                                         family='freeserif', color='#000099', )
            ax1 = fig.add_subplot(2, 2, j + 1)
            ax1.xaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.1, zorder=0)
            ax1.yaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.1, zorder=1)
            #print counter
            ax1.set_title(str(int(earthquakes[counter, 1])) + '/' + str(int(earthquakes[counter, 2])) + '/' +
                          str(int(earthquakes[counter, 3])))
            ax1.loglog(freqs, momentrate[:, counter], lw=2, ls='-', color='#FF8000')
            ax1.loglog(freqs, momentratepstd[:, counter], lw=2, ls='--',
                       zorder=2, color='#193300')
            ax1.loglog(freqs, momentratenstd[:, counter], lw=2, ls='--',
                       zorder=2, color='#193300')
            fit = nppoly.polyfit(np.log(freqs), np.log(momentrate[:, counter]), 1)
            ax1.loglog(freqs, np.exp(nppoly.polyval(np.log(freqs), fit)),
               alpha=0.5, color='#994c00', linewidth=2,
               label=r'Fitted line  $Q_{s}=%0.2f \times f^{%0.1f}$' % (50, fit[1]), zorder=3)
            counter += 1
            # ax1.set_ylim([10**20, 10**26])
            ax1.set_xlim([0, 100])
        plt.tight_layout()
        fig.savefig(outputpath + 'Source-effects/Source-plots'
                    + str(i) + '.pdf', format='pdf', dpi=200)
        plt.close(fig)


def gridsearch(inputpath, outputpath, sigma, gamma, magnitude, bs=3.5, n=20, ps=2.8, alphas=6):
    earthquakes = np.genfromtxt(inputpath + 'earth.csv',
                                delimiter=',')
    freqs = np.genfromtxt(outputpath + 'Interpolating/freqs.txt',
                          delimiter=',')
    momentrates = np.genfromtxt(outputpath +
                                'Source-effects/momentrates.txt', delimiter=',')
    magnitudes = np.genfromtxt(outputpath +
                               'Source-effects/magnitudesextracted.txt', delimiter=',')

# --------------   Introducing variables --------------------
    momentratescalculated = np.zeros([(np.shape(momentrates)[0]), (np.shape(momentrates))[1]])
    gammanumberofsamples = gamma[2]
    sigmaanumberofsamples = sigma[2]
    magnitudenumberofsamples = magnitude[1]
    dsigma = np.linspace(sigma[0], sigma[1], sigmaanumberofsamples)
    dgamma = np.linspace(gamma[0], gamma[1], gammanumberofsamples)
    magnitudessamples = np.linspace(-magnitude[0], magnitude[0], magnitudenumberofsamples)
    bestanswer = np.zeros([len(magnitudes), 5], dtype=float)
    bestanswer2 = np.zeros([len(magnitudes), 5], dtype=float)
    forplot = np.zeros([len(magnitudes), 4], dtype=float)
    momentrateaverage = np.zeros(len(magnitudes))
    momentrateaverage = np.mean(momentrates, axis=0)  # Average of moment rate in each record
    swaveenergy = np.zeros(len(magnitudes))
    energyconstant = ((15 * np.pi * 1000 * ps * ((alphas * 1000) ** 5)) ** (-1)) + ((15 * np.pi * ps * 1000 *((bs *1000) ** 5)) ** (-1)) # Constant for S-wave ebergy calculation
    swaveenergy = energyconstant * np.asarray([np.average([(2 * np.pi * freqs[i] * momentrates[i, z]) ** 2 for i in range(n)])
                                     for z in range(len(magnitudes))])  # S-wave energy calculations
    bestanswereachfreq = np.zeros([len(magnitudes), 6], dtype=float)
    earthquakeprint = np.zeros([len(magnitudes)], dtype='S10, float, float, float, float, float, float, float, float')
    bestanswer[:, 0] = 10000000
    bestanswer2[:, 0] = 10000000
    #print magnitudessamples
# ------------------------- Grid search process ---------------------------
    fc = 0
    f = open(outputpath +
             'Source-effects/grid-search/sigmagammacomparisons.txt', 'w')
    gg = open(outputpath +
              'Source-effects/grid-search/Comparison-for-each-frequency.txt', 'w')
    for z in range(len(magnitudes)):
        difference = 100000000000
        print >> gg, '\n\n For earthquake number %d and %0.2fM magnitude :  /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/' \
                     % (z + 1, magnitudes[z])
        bestanswer[z, 0] = 10000000.0

        for p in range(magnitudenumberofsamples):  # forloop for magnitudes calculations
            momentscalculated = 10 ** (1.5 * (magnitudes[z] + magnitudessamples[p] + 10.7))

            for k in range(sigmaanumberofsamples):  # forloop for sigma calculations
                fc = 4.9 * (10 ** 6) * bs * ((dsigma[k] / momentscalculated) ** (1.0 / 3.0))

                for j in range(gammanumberofsamples):  # forloop for gamma calculations
                    momentratescalculated[:, z] = np.log(momentscalculated / (1 + (((freqs[:]) / fc) ** dgamma[j])))
                    momentratelog = np.log(momentrates)
                    difference = (np.mean((momentratelog[: ,z] - momentratescalculated[:, z]) ** 2))

                    if difference < bestanswer[z, 0]:
                        bestanswer[z, 0] = (difference)
                        bestanswer[z, 1] = dgamma[j]
                        bestanswer[z, 2] = dsigma[k]
                        bestanswer[z, 3] = magnitudes[z] + magnitudessamples[p]
                        bestanswer[z, 4] = np.average(np.exp(momentratescalculated[:, z]))
                        fcc = fc

        print >> gg, 'For %0.2f Hz, best answer is : \n ' \
                     'gamma = %0.2f, Sigma = %0.2f, magnitude = %0.2f, moment = %0.2e, fc=%0.2f Difference = %0.12f   -\n' \
                     % (bestanswereachfreq[i, 0], bestanswereachfreq[i, 2], bestanswereachfreq[i, 3],
                        bestanswereachfreq[i, 4], bestanswereachfreq[i, 5], fc, bestanswereachfreq[i, 1])
                    # print np.shape(bestanswereachfreq)
                    # print np.hstack((bestanswereachfreq, bestanswer))
        forplot[z, 0] = (2.34 * bs) / (2 * np.pi * fcc)    # Radius
        forplot[z, 1] = fcc  #corner frequency
        forplot[z, 2] = swaveenergy[z]
        forplot[z, 3] = momentrateaverage[z]
        earthquakeprint[z] = (str(int(earthquakes[z, 1])) + '/' +
                              str(int(earthquakes[z, 2])) + '/' + str(int(earthquakes[z, 3])),
                              bestanswer[z, 0], bestanswer[z, 1], bestanswer[z, 2], bestanswer[z, 3],
                              magnitudes[z], fcc, forplot[z,0], bestanswer[z, 4])

    f.close()
    gg.close()
    textheader = "Earthquake date, Difference, Dgamma, Dsigma,New magnitude, Old magnitude, fc, radius, New momentrate"
    np.savetxt(outputpath +
               'Source-effects/grid-search/gridsearch-results.txt', earthquakeprint,
               header=textheader, fmt='%-10s, %10.2e, %5.2f, %10.2f, %5.2f, %5.2f, %5.2f, %5.2f, %10.2e',
               delimiter=',')

def gridsearchplots(inputpath, outputpath, eqcount, stcount, n=20 ):

    gridresults = np.genfromtxt(outputpath + 'Source-effects/grid-search/gridsearch-results.txt',
                                delimiter=',')
    earthquakes = np.genfromtxt(inputpath + 'earth.csv',
                                delimiter=',')
    freqs = np.genfromtxt(outputpath + 'Interpolating/freqs.txt',
                          delimiter=',')
    momentrate = np.genfromtxt(outputpath +
                                'Source-effects/momentrates.txt', delimiter=',')
    momentratepstd = np.genfromtxt(outputpath +
                                'Source-effects/momentratespstd.txt', delimiter=',')
    momentratenstd = np.genfromtxt(outputpath +
                                'Source-effects/momentratesnstd.txt', delimiter=',')
    magnitudes = np.genfromtxt(outputpath +
                               'Source-effects/magnitudesextracted.txt', delimiter=',')
    momentratecalculated = np.zeros([n, eqcount])
    #print gridresults
    #print freqs[2]
    momentratecalculated = [(10 ** (1.5 * (gridresults[z, 4] + 10.7)) / (1 + (f / gridresults[z, 7]) ** gridresults[z, 2])) for z in range(int(eqcount)) for f in np.linspace(0.4, 15, 1000)]
    #print np.shape(momentratecalculated)
    counter = 0
    for i in range(int(eqcount) / 4 + 1):
        if counter == eqcount:
            break
        fig = plt.figure(i)
        for j in range(4):
            if counter == eqcount:
                break
            ax1 = fig.add_subplot(2, 2, j + 1)
            if j == 0 or j == 2:
                ax1.yaxis.set_label_text(r'Moment rate spectrum $(dyne \times Cm)$', size=12,
                                         family='freeserif', color='#000099', )
            ax1 = fig.add_subplot(2, 2, j + 1)
            ax1.xaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.1, zorder=0)
            ax1.yaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.1, zorder=1)
            ax1.set_title(str(int(earthquakes[counter, 1])) + '/' + str(int(earthquakes[counter, 2])) + '/' +
                          str(int(earthquakes[counter, 3])))
            ax1.loglog(freqs, momentrate[:, counter], lw=2, ls='-', color='#FF8000')
            ax1.loglog(freqs, momentratepstd[:, counter], lw=2, ls='--',
                       zorder=2, color='#193300')
            ax1.loglog(freqs, momentratenstd[:, counter], lw=2, ls='--',
                       zorder=2, color='#193300')
            fit = nppoly.polyfit(np.log(freqs), np.log(momentrate[:, counter]), 1)
            ax1.loglog(freqs, np.exp(nppoly.polyval(np.log(freqs), fit)),
               alpha=0.5, color='#994c00', linewidth=2)
            ax1.loglog(np.linspace(0.4, 15, 1000), momentratecalculated[counter*1000: counter*1000 + 1000],
               alpha=0.5, color='#994c00', linewidth=2)
            counter += 1
            # ax1.set_ylim([10**20, 10**26])
            ax1.set_xlim([0, 100])
        plt.tight_layout()
        fig.savefig(outputpath + 'Source-effects/grid-search/source-plot'
                    + str(i) + '.pdf', format='pdf', dpi=200)
        plt.close(fig)

def extrasourceplots(inputpath, outputpath, eqcount, stcount, n=20, ps=2.8, bs=3.5, alphas=5.98):

    gridresults = np.genfromtxt(outputpath + 'Source-effects/grid-search/gridsearch-results.txt',
                                delimiter=',')
    earthquakes = np.genfromtxt(inputpath + 'earth.csv',
                                delimiter=',')
    freqs = np.genfromtxt(outputpath + 'Interpolating/freqs.txt',
                          delimiter=',')
    momentrates = np.genfromtxt(outputpath +
                                'Source-effects/momentrates.txt', delimiter=',')
    momentratepstd = np.genfromtxt(outputpath +
                                'Source-effects/momentratespstd.txt', delimiter=',')
    momentratenstd = np.genfromtxt(outputpath +
                                'Source-effects/momentratesnstd.txt', delimiter=',')
    magnitudes = np.genfromtxt(outputpath +
                               'Source-effects/magnitudesextracted.txt', delimiter=',')
    momentratecalculated = np.zeros([n, eqcount])

    momentrateaverage = np.zeros(len(magnitudes))
    momentrateaverage = np.mean(momentrates, axis=0)  # Average of moment rate in each record
    momentscalculatedfromgrid = 10 ** (1.5 * (gridresults[:, 4] + 10.7))

    swaveenergy = np.zeros([n, len(magnitudes)])
    swaveenergy2 = np.zeros(len(magnitudes))
    energyconstant = ((15 * np.pi * 1000 * ps * ((alphas * 1000) ** 5)) ** (-1)) + ((10 * np.pi * ps * 1000 *((bs *1000) ** 5)) ** (-1)) # Constant for S-wave ebergy calculation

    for z in range (len(magnitudes)) :
        swaveenergy[:, z] = (2 * np.pi * freqs[:] * momentrates[:, z] * 10 ** (-7)) ** 2
        swaveenergy2[z] = energyconstant * np.trapz(swaveenergy[:, z], freqs) * 10 ** (7)

    radius = gridresults[:, 7]
    radiussample = np.linspace(0.1, 100, 200)
    fc = gridresults[:, 6]
# ------------------------- Radius - Moment Plot --------------------------------

    figoor = plt.figure(122)
    ax2 = figoor.add_subplot(111)
    ax2.scatter(radius, momentscalculatedfromgrid, s=100,
                color='deepskyblue', edgecolor='black', alpha=0.5, zorder=3,)
    ax2.loglog(radiussample, 1 * 16 * ((radiussample) ** 3)* (10 ** 21) * (1.0/7.0),
               color = 'peachpuff', linewidth = 3, zorder=1, label='1 bar')
    ax2.loglog(radiussample, 10 * 16 * ((radiussample) ** 3)* (10 ** 21) * (1.0/7.0),
               color = 'lightsalmon', linewidth = 3, label='10 bar')
    ax2.loglog(radiussample, 100 * 16 * ((radiussample) ** 3)* (10 ** 21) * (1.0/7.0),
               color = 'coral', linewidth = 3, label='100 bar')
    ax2.loglog(radiussample, 1000 * 16 * ((radiussample) ** 3)* (10 ** 21) * (1.0/7.0),
               color = 'red', linewidth = 3, label='1000 bar')
    ax2.set_xlim([0, 20])
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.xaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.3, zorder=0)
    ax2.yaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.3, zorder=1)
    ax2.xaxis.set_label_text('Radius', size=18, family='freeserif', color='#000099')
    ax2.yaxis.set_label_text('Seismic moment (dyne.cm)', size=18, family='freeserif',
                             color='#000099')
    ax2.legend(loc='upper left', frameon=False, fancybox=True, shadow=True,
               prop={'size': 18, 'family': 'freeserif'})
    figoor.savefig(outputpath + 'Source-effects/Extra-plots/'
                    + 'radius-moment' + '.pdf', format='pdf', dpi=200)


# -------------------------  Moment - S-wave energy Plot --------------------------------

    fig = plt.figure(214)
    ax = fig.add_subplot(111)
    ax.scatter(momentscalculatedfromgrid, swaveenergy2, s=100,
               color='deepskyblue', edgecolor='black', alpha=0.5, zorder=3,)
    ax.xaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.1, zorder=0)
    ax.yaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.1, zorder=1)
    fit = nppoly.polyfit(np.log(momentscalculatedfromgrid), np.log(swaveenergy2), 1)
    fit2 = nppoly.polyfit(np.log(momentscalculatedfromgrid), np.log(swaveenergy2/momentscalculatedfromgrid), 0)

    ax.loglog(momentscalculatedfromgrid, np.exp(nppoly.polyval(np.log(momentscalculatedfromgrid), fit)),
               alpha=0.5, color='#994c00', linewidth=2,
               label='Fitted line', zorder=3)
    fig.suptitle(r'$\frac{E_{s}}{M} = %0.2e$' % np.exp(fit2))
    ax.xaxis.set_label_text('Seismic Moment (dyne.Cm)', size=16, family='freeserif', color='#000099')
    ax.yaxis.set_label_text('S-wave Energy (dyne.Cm)', size=16, family='freeserif',
                             color='#000099')
    # ax.set_xlim([0.1, 100])
    ax.set_xscale('log')
    ax.set_yscale('log')
    fig.savefig(outputpath + 'Source-effects/Extra-plots/'
                   + 'Es-moment' + '.pdf', format='pdf', dpi=200)


# ------------------------- Moment - fc ---------------------------------

    fig = plt.figure(224)
    ax = fig.add_subplot(111)
    ax.scatter(fc, momentscalculatedfromgrid, s=100,
               color='deepskyblue', edgecolor='black', alpha=0.5, zorder=3,)
    ax.xaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.1, zorder=0)
    ax.yaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.1, zorder=1)
    fit = nppoly.polyfit(np.log(fc), (np.log(momentscalculatedfromgrid) + 3 * np.log(fc)), 0)
    ax.loglog(fc, np.exp(-3 * np.log(fc) + fit),
               alpha=0.5, color='#994c00', linewidth=2,
               label='Fitted line', zorder=3)

    regionaldsigma = np.exp(fit)/(4.9 * 10 ** 6 * bs) ** 3
    fig.suptitle(r'$\delta \sigma = %0.2f$ bar' % regionaldsigma)
    ax.xaxis.set_label_text('Corner Frequecy (Hz)', size=18, family='freeserif', color='#000099')
    ax.yaxis.set_label_text('Seismic moment (dyne.Cm)', size=18, family='freeserif',
                             color='#000099')
    ax.set_xlim([0.1, 100])
    ax.set_xscale('log')
    ax.set_yscale('log')
    fig.savefig(outputpath + 'Source-effects/Extra-plots/'
                   + 'fc-moment' + '.pdf', format='pdf', dpi=200)


# ------------------------  M - Distance plot -----------------------------

    spec = np.genfromtxt(inputpath +
                         '/spectrums.csv', delimiter=',')

    fig = plt.figure(334)
    ax = fig.add_subplot(111)
    dismag = [gridresults[spec[0, i]-1, 4] for i in range(np.shape(spec)[1])]
    ax.scatter(spec[2, :], dismag, s=100,
               color='deepskyblue', edgecolor='black', alpha=0.5, zorder=3,)
    ax.xaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.1, zorder=0)
    ax.yaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.1, zorder=1)

    fig.suptitle('Distance - Magnitude')
    ax.xaxis.set_label_text('Hypo Distance (Km)', size=18, family='freeserif', color='#000099')
    ax.yaxis.set_label_text('Magnitude (Mw)', size=18, family='freeserif',
                             color='#000099')

    fig.savefig(outputpath + 'Source-effects/Extra-plots/'
                   + 'M-distance' + '.pdf', format='pdf', dpi=200)