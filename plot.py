import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch

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

grdata = pd.read_fwf('gr.dat',skiprows=1, widths=[8,12],header=None)
grdata.columns = ['r','gr']

plt.figure(figsize=(6,4))
plt.plot(grdata.r, grdata.gr, lw=3,label='a-Si model')
plt.xlabel('r ($\AA$)')
plt.ylabel('gr($r$)')
plt.xlim([1,20])
plt.ylim([0,6])
plt.title('Pair correlatin function')
plt.legend()
plt.savefig('gr.png', dpi=100, transparent=False, bbox_inches='tight', pad_inches=0.1)
#plt.show()

skdata = pd.read_fwf('Sk.dat',skiprows=0, widths=[8,12],header=None)
skdata.columns = ['k','sk']

plt.figure(figsize=(6,4))
plt.plot(skdata.k, skdata.sk, 'r--o', ms=5,lw=1,label='a-Si model',alpha=.2)
plt.xlabel('k ($\AA^{-1}$)')
plt.ylabel('S($k$)')
plt.xlim([.5,25])
plt.ylim([0,2.2])
plt.title('Structure Factor')
plt.legend()
plt.savefig('Sk.png', dpi=100, transparent=False, bbox_inches='tight', pad_inches=0.1)
#plt.show()

angledata = pd.read_fwf('bad.dat',skiprows=1, widths=[12,12,12],header=None)
angledata.columns = ['ang','dist','ndist']

plt.figure(figsize=(6,4))

xx=angledata.ang.to_numpy()
yy=angledata.ndist.copy().to_numpy()

path = Path(np.array([xx,yy]).transpose())
patch = PathPatch(path, facecolor='none')
plt.gca().add_patch(patch)

plt.plot(xx,yy,lw=3)
plt.xlabel('$\Theta$')
plt.ylabel('Distribution')
plt.xlim([60,160])
plt.ylim([0,.3])
plt.title('Bond-Angle Distribution')
plt.legend(['aSi model'])

im=plt.imshow(xx.reshape(yy.size,1),  cmap=plt.cm.Blues,interpolation="bicubic",
                origin='lower',extent=[60,170,0.0,.3],aspect="auto", clip_path=patch, clip_on=True)
plt.savefig('bad.png', dpi=100, transparent=False, bbox_inches='tight', pad_inches=0.1)
#plt.show()

