import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

plt.rcParams['font.family'] = 'serif'
plt.rcParams["mathtext.fontset"] = "dejavuserif"

source_name = '3C266'
wav = '450'
n = 250

snr = fits.open('gs20201102_78_850_snr_crop.fit')[0].data[0].flatten()
# jks = fits.open(f'reduced-jk/mosaic/{source_name}_{wav}_jk_crop_mf_cal_snr.fits')[0].data[0].flatten()
snr = snr[~np.isnan(snr)]
# jks = jks[~np.isnan(jks)]

rang = np.hstack([snr])
rang = np.nanmin(rang),np.nanmax(rang)

fig,ax = plt.subplots(dpi=120)

s = ax.hist(snr,range=rang,bins=n,histtype='step',label='Original')
# ax.hist(jks,range=rang,bins=n,histtype='step',label='Jackknife')
ax.set_xlabel('SNR')
ax.set_xlim(-max(s[1])-1,max(s[1])+1)
ax.set_yscale('log')
ax.set_ylabel('Counts')
ax.set_title(rf'{source_name} ${wav}\mu m$')
ax.legend()

#plt.savefig(f'{source_name}_{wav}_flux dist.png')
plt.show()
