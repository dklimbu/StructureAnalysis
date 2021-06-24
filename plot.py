import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

grdata = pd.read_fwf('gr.dat',skiprows=0, widths=[8,12],header=None)
grdata.columns = ['r','gr']


parameters = {'font.size': 15,
                    'axes.labelsize': 20,
                    'xtick.labelsize': 18,
                    'ytick.labelsize': 18,
                    'legend.fontsize': 15,
                    'axes.linewidth': 2,
                    'lines.linewidth': 2,
                    'xtick.major.width': 3,
                    'ytick.major.width': 3,}
plt.rcParams.update(parameters)

plt.figure(figsize=(6,4))
plt.plot(grdata.r, grdata.gr, lw=3)
plt.xlabel('r ($\AA$)')
plt.ylabel('gr($r$)')
plt.xlim([1,20])
plt.ylim([0,6])
plt.title('Pair correlatin function')
plt.savefig('gr.png', dpi=100, transparent=False, bbox_inches='tight', pad_inches=0.1)

#plt.show()


angledata = pd.read_fwf('bad.dat',skiprows=1, widths=[12,15],header=None)
angledata.columns = ['ang','dist']


plt.figure(figsize=(6,4))
plt.plot(angledata.ang, angledata.dist/5000, lw=3)
plt.xlabel('$\Theta$')
plt.ylabel('Distribution')
plt.xlim([60,160])
plt.ylim([0,.3])
plt.title('Bond-Angle Distribution')
plt.savefig('bad.png', dpi=100, transparent=False, bbox_inches='tight', pad_inches=0.1)

#plt.show()

