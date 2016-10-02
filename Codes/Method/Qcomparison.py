import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 1.75
mpl.rcParams['font.size'] = 12
mpl.rcParams['font.family'] = 'serif'
import numpy as np
import matplotlib.pyplot as plt


def QComp(f1, f2, n):
    a1 = (np.log(f2) - np.log(f1)) / (n - 1)
    w1 = np.arange(0, n, dtype=float)
    w1[0] = f1
    for i in range(1, int(n)):
        w1[i] = np.exp(np.log(w1[0])+(i*a1))
    #  -------------------- Q models-------------------
    weak = [183 * f ** 0.61 for f in w1]
    strong = [186 * f ** 0.66 for f in w1]
    ghassemi = [90 * f ** 0.74 for f in w1]
    mousavi = [153 * f ** 0.88 for f in w1]
    Hassani = [151 * f ** 0.75 for f in w1]

    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)
    ax1.loglog(w1, weak, c='b', lw='2', marker='D', label = 'Weak-motion Data used')
    ax1.loglog(w1, strong, c='black', lw='3', ls='-', label = 'This study')
    ax1.loglog(w1, ghassemi, c='orange', lw='2', ls=':', label = 'Previous studies-"Ghassemi et al"')
    ax1.loglog(w1, mousavi, c='purple', lw='2', marker='*', label = 'Previous studies-"Mousavi  et al"')
    ax1.loglog(w1, Hassani, c='brown', lw='2', marker='o', label = 'Previous studies-"Hassani et al"')


    ax1.xaxis.set_label_text('Frequency(Hz)', size=12, family='freeserif', color='#000099')
    ax1.yaxis.set_label_text('Q-factor(Qs)', size=12, family='freeserif', color='#000099')
    ax1.xaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.5, zorder=0)
    ax1.yaxis.grid(which='both', ls='-', lw=0.3, color='#c0c0c0', alpha=0.5, zorder=1)
    ax1.legend(loc='upper left', frameon=False, fancybox=True, shadow=True,
               prop={'size': 10, 'family': 'freeserif'})
    ax1.set_xlim([0.4, 30])
    ax1.set_ylim([0, 1000])
    plt.show()
    fig.savefig('/home/babak/PycharmProjects/Simulation/GeneralInversion/Results/Path-effect/Qcomparison.pdf',
                format='pdf', dpi=1000)


QComp(0.4, 15, 20)