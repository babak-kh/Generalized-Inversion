import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['grid.color'] = '#6666FF'
import matplotlib.pyplot as plt
import numpy as np
import random

fig = plt.figure(1, figsize=(10,5), facecolor='white', edgecolor='red')
fig.suptitle('This is the end', fontsize=14, fontweight='bold',fontname='freeserif', color='red')
#fig.set_suptitle('NIIICE')


ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=2, rowspan=1)
ax1.xaxis.set_label_text('salam')
ax1.yaxis.set_label_text('chakerim', size=10, family='Rekha', style='italic')
ax1.tick_params(axis='both', labelcolor='red', gridOn='True', size =2, width=5)
ax4 = plt.subplot2grid((3, 3), (0, 2), colspan=1, rowspan=1, sharey=ax1, sharex=ax1)
ax2 = plt.subplot2grid((3, 3), (1, 0), colspan=3, rowspan=1)
ax2.tick_params(axis='both', labelcolor='red', gridOn='True', size =2, width=5)


ax3 = plt.subplot2grid((3, 3), (2, 0), colspan=3, rowspan=1)
ax3.xaxis.set_label_text('XaXiS', fontsize='10', color='#000099', family='freeserif')
ax3.yaxis.set_label_text('YaXiS $salam \pi \omega \cos 2$ ', fontsize='10', color='red', family='freeserif')

ax3.xaxis.grid(which='both',lw=1, ls='--')









ax1.plot(np.arange(100), linewidth=1, color='red', linestyle='--')
ax1.grid(True)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax2.plot(np.arange(200))
ax3.plot(np.arange(500))

#plt.tight_layout()
plt.show()