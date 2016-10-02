import numpy as np
import matplotlib.pyplot as plt
import numpy.polynomial.polynomial as nppoly


outputpath = '/home/babak/PycharmProjects/generalwithinterface/GeneralInversion/Results/'
n = 20

result = np.genfromtxt(outputpath + 'svd/answer.txt',
                       delimiter=',')
covd = np.genfromtxt(outputpath + 'svd/Covariance/covd.txt',
                     defaultfmt='.4e', delimiter=',')
freqs = np.genfromtxt(outputpath + 'Interpolating/freqs.txt',
                      delimiter=',')
error = np.zeros(n)
qfactor2 = result[:n]
for j in range(n):
    s = [np.random.normal(result[j], np.sqrt(covd[j,j])) for x in range(1000)]
    s = np.asarray(s)
    s = 1/s
    error[j] = np.std(s)
qfactor = 1.0 / qfactor2

np.savetxt(outputpath + 'Path-effect/q.txt', qfactor,
           fmt='%0.5f', delimiter=',')
np.savetxt(outputpath + 'Path-effect/covdiagonal.txt', np.sqrt(np.diag(covd)[:]),
           fmt='%0.5f', delimiter=',')
fit = nppoly.polyfit(np.log(freqs), np.log(qfactor), 1)

aa = np.exp(fit[0])

# ----------------------------- plotting -------------------

fig = plt.figure(1)
fig.suptitle('Q-factor', fontsize=14, fontweight='bold', family='freeserif')
ax1 = fig.add_subplot(111)
ax1.xaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.1, zorder=0)
ax1.yaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.1, zorder=1)
ax1.fill_between(freqs, qfactor+error+0.5, qfactor-error-0.5, facecolor = 'lightblue', alpha = 0.5 )
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

ax1.set_xlim([0.2, 20])
ax1.set_ylim([90, 1500])
ax1.legend(loc='upper left', frameon=False, fancybox=True, shadow=True,
           prop={'size': 14, 'family': 'freeserif'})
fig.savefig(outputpath + 'Path-effect/Quality-factor.pdf', format='pdf', dpi=1000)
plt.close(fig)
