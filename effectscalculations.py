import matplotlib as mpl

mpl.rcParams['lines.linewidth'] = 1.75
mpl.rcParams['font.size'] = 12
mpl.rcParams['font.family'] = 'serif'
import numpy as np
from GeneralInversion.Codes import GeneralizedInversion
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as nppoly
import math as mt
import xlrd


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
            efromsvd[j, i] = result[n + (i * n) + j]
            momentrate[j, i] = np.exp(result[n + (i * n) + j]) / (4 * c * (np.pi ** 2) * (freqs[j] ** 2))
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
               momentrate, fmt='%0.5f', delimiter=',')
    np.savetxt(outputpath + 'Source-effects/magnitudesextracted.txt',
               moment, fmt='%0.5f', delimiter=',')
    counter = 0
    for i in range(eqcount / 4 + 1):
        fig = plt.figure(i)
        for j in range(4):
            if (j + 1) * (i + 1) > eqcount:
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
            counter += 1
            # ax1.set_ylim([10**20, 10**26])
            ax1.set_xlim([0, 100])
        plt.tight_layout()
        fig.savefig(outputpath + 'Source-effects/Source-plots'
                    + str(i) + '.pdf', format='pdf', dpi=200)
        plt.close(fig)


def gridsearch(inputpath, outputpath, sigma, gamma, magnitude, bs=3.5, n=20):
    earthquakes = np.genfromtxt(inputpath + 'earth.csv',
                                delimiter=',')
    freqs = np.genfromtxt(outputpath + 'Interpolating/freqs.txt',
                          delimiter=',')
    momentrates = np.genfromtxt(outputpath +
                                'Source-effects/momentrates.txt', delimiter=',')
    magnitudes = np.genfromtxt(outputpath +
                               'Source-effects/magnitudesextracted.txt', delimiter=',')
    # print gamma
    # print magnitudes
    gammanumberofsamples = gamma[2]
    sigmaanumberofsamples = sigma[2]
    magnitudenumberofsamples = magnitude[1]
    dsigma = np.linspace(sigma[0], sigma[1], sigmaanumberofsamples)
    dgamma = np.linspace(gamma[0], gamma[1], gammanumberofsamples)
    magnitudessamples = np.linspace(-magnitude[0], magnitude[0], magnitudenumberofsamples)
    # ----------  test  --------------------------------
    # momentscalculated = 10 ** (1.5 * (magnitudes[0] + 0.2 + 10.7))
    # fc = 4.9 * 10**6 * bs * ((200 / momentscalculated) ** (1.0/3.0))
    # momentratescalculated = momentscalculated / (1+(((freqs[0]) / fc) ** 2))
    # # momentratescalculated = np.log(momentratescalculated)
    # difference = np.log(momentrates[0, 0]) - np.log(momentratescalculated)
    # print momentrates[0,0]
    # print momentratescalculated
    # print difference
    # ------------------------  test  --------------------
    bestanswer = np.zeros([len(magnitudes), 5], dtype=float)
    bestanswereachfreq = np.zeros([len(magnitudes), 6], dtype=float)
    earthquakeprint = np.zeros([len(magnitudes)], dtype='S10, float, float, float, float, float, float, float')
    bestanswer[:, 0] = 10000000
    fc = 0
    f = open(outputpath +
             'Source-effects/grid-search/sigmagammacomparisons.txt', 'w')
    gg = open(outputpath +
              'Source-effects/grid-search/Comparison-for-each-frequency.txt', 'w')
    for z in range(len(magnitudes)):
        print >> f, 'For earthquake number %d with %0.2fM magnitude ---------------\n\n' \
                    % (z + 1, magnitudes[z])
        print >> gg, 'For earthquake number %d and %0.2fM magnitude :  /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/' \
                     % (z + 1, magnitudes[z])
        bestanswer[z, 0] = 10000000.0

        for i in range(n):  # forloop for comparison between the moments
            for j in range(gammanumberofsamples):  # forloop for gamma calculations
                for k in range(sigmaanumberofsamples):  # forloop for sigma calculations
                    for p in range(magnitudenumberofsamples):  # forloop for magnitudes calculations

                        momentscalculated = 10 ** (1.5 * (magnitudes[z] + magnitudessamples[p] + 10.7))
                        fc = 4.9 * 10 ** 6 * bs * ((dsigma[k] / momentscalculated) ** (1.0 / 3.0))
                        momentratescalculated = momentscalculated / (1 + (((freqs[i]) / fc) ** dgamma[j]))
                        difference = np.mean(np.log(momentrates[i, z]) - np.log(momentratescalculated)) ** 2
                        print >> f, 'for gamma=%0.2f and sigma=%0.2f and magnitude=%0.2f in %0.2f Hz, difference is %0.5f' \
                                    % (
                                    dgamma[j], dsigma[k], (magnitudes[z] + magnitudessamples[p]), freqs[i], difference)
                        if abs(difference) < bestanswer[z, 0]:
                            bestanswer[z, 0] = abs(difference)
                            bestanswer[z, 1] = dgamma[j]
                            bestanswer[z, 2] = dsigma[k]
                            bestanswer[z, 3] = magnitudes[z] + magnitudessamples[p]
                            bestanswer[z, 4] = momentratescalculated
            bestanswereachfreq[i, :] = np.hstack((freqs[i], bestanswer[z, :]))
            print >> gg, 'For %0.2f Hz, best answer is : \n ' \
                         'gamma = %0.2f, Sigma = %0.2f, magnitude = %0.2f, moment = %0.2e, Difference = %0.12f   -\n' \
                         % (bestanswereachfreq[i, 0], bestanswereachfreq[i, 2], bestanswereachfreq[i, 3],
                            bestanswereachfreq[i, 4], bestanswereachfreq[i, 5], bestanswereachfreq[i, 1])
            # print np.shape(bestanswereachfreq)
            # print np.hstack((bestanswereachfreq, bestanswer))

        earthquakeprint[z] = (str(int(earthquakes[z, 1])) + '/' +
                              str(int(earthquakes[z, 2])) + '/' + str(int(earthquakes[z, 3])),
                              bestanswer[z, 0], bestanswer[z, 1], bestanswer[z, 2], bestanswer[z, 3],
                              magnitudes[z], fc, bestanswer[z, 4])
    f.close()
    textheader = "Earthquake date, Difference, Dgamma, Dsigma,Old magnitude, New magnitude, fc, New momentrate"
    np.savetxt(outputpath +
               'Source-effects/grid-search/gridsearch-results.txt', earthquakeprint,
               header=textheader, fmt='%-10s, %10.2e, %5.2f, %10.2f, %5.2f, %5.2f, %5.2f, %10.2e',
               delimiter=',')


def gridsearchplots(inputpath, outputpath):
    sourcemomentrates = np.genfromtxt(outputpath + 'Source-effects/momentrates.txt',
                                      delimiter=',')
    gridresults = np.genfromtxt(outputpath + 'Source-effects/grid-search/gridsearch-results.txt',
                                delimiter=',')
    freqs = np.genfromtxt(outputpath + 'Interpolating/freqs.txt',
                          delimiter=',')
    # for i in range(np.shape(sourcemomentrates)[1]):
    #     fig = plt.figure(1)
    #     ax1 = fig.add_subplot(111)
    #     ax1.loglog(freqs, sourcemomentrates[:, i])
    #     for j in range(len(freqs)):
    #         a[j] = 4.9 * 10 ** 6 * 3.5 *
    #     ax1.loglog(freqs, )

