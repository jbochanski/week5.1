# -*- coding: utf-8 -*-
#creates templates of M dwarfs and M giants and overplots our science target
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
import asciitable
from matplotlib.ticker import MultipleLocator

#mpl.rcParams['lines.linewidth']= 2.0
#mpl.rcParams['patch.edgecolor']=   'white'
#mpl.rcParams['axes.color_cycle']=   ['377EB8', 'E41A1C', '4DAF4A', '984EA3', 'FF7F00','#FFFF33','#A65628', '#F781BF', '#999999']
#mpl.rcParams['image.cmap']=  'gray'

#mpl.rcParams['xtick.major.size'] = 16
#mpl.rcParams['xtick.major.width'] = 2
#mpl.rcParams['xtick.minor.size'] = 9
#mpl.rcParams['xtick.minor.width'] = 2
#mpl.rcParams['ytick.major.size'] = 16
##mpl.rcParams['ytick.major.width'] = 2
#mpl.rcParams['ytick.minor.size'] = 9
#mpl.rcParams['ytick.minor.width'] = 2
#mpl.rc('xtick', labelsize=25) 
#mpl.rc('ytick', labelsize=25) 
#mpl.rcParams['xtick.major.pad'] = 8
#mpl.rcParams['ytick.major.pad'] = 8

smooth_beta = 4

class Spectrum(object):
    def __init__(self, wave=None, flux=None):
        self.wave = wave
        self.flux = flux
        
              
def read_mmt_spec(filename):
    d1 = fits.open(filename)
    fluxd = d1[0].data
    c0 = d1[0].header['CRVAL1']
    c1 = d1[0].header['CDELT1']
    npix = d1[0].header['naxis1']
    waved = 10.**(c0 + c1 * np.arange(npix))
    z = np.polyfit(waved, fluxd, 4)
    f = np.poly1d(z)
    fluxd = fluxd/f(waved)
    #median(fluxd[g1])
    data = Spectrum(waved,fluxd)
    return data
    
def read_lbt_spec(filename):
    d= asciitable.read(filename, names=('wave','flux','err'))
    waved = np.array(d['wave'])
    fluxd = np.array(d['flux'])
    g1 = ((waved > 7000) & (waved < 9000))
    waved = waved[g1]
    fluxd = fluxd[g1]
    z = np.polyfit(waved, fluxd, 4)
    f = np.poly1d(z)
    fluxd = fluxd/f(waved)
    #median(fluxd[g1])
    data = Spectrum(waved,fluxd)
    return data

def smooth(x,beta):
    #""" kaiser window smoothing """
    window_len=11
 # extending the data at beginning and at the end
 # to apply the window at the borders
    s = np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    w = np.kaiser(window_len,beta)
    y = np.convolve(w/w.sum(),s,mode='valid')
    return y[5:len(y)-5]
    
    

sci_spec_file = '/Users/jjb/Dropbox/ukidss_m_giants/data/mmt/final/Mg3.90_1.93_f.fits'
sci_spec_file2 = '/Users/jjb/Dropbox/ukidss_m_giants/data/mmt/final/Mg116.07_25.54_f.fits'

sci_spec = read_mmt_spec(sci_spec_file)
sci_spec2 = read_mmt_spec(sci_spec_file2)


giant_spec_file = '/Users/jjb/Dropbox/ukidss_m_giants/data/mmt/final/HIP13756_f.fits'
#giant_spec_file = '/Users/jjb/Dropbox/ukidss_m_giants/data/mmt/final/HD5820_f.fits'  #using this even though its not the best chi-sq match

giant_spec = read_mmt_spec(giant_spec_file)
dwarf_spec_file = '/Users/jjb/Dropbox/ukidss_m_giants/data/mmt/final/GXAnd_f.fits'

#dwarf_spec_file = '/Users/jjb/Dropbox/ukidss_m_giants/data/mmt/final/HIP9291_f.fits'
dwarf_spec = read_mmt_spec(dwarf_spec_file)



sci_spec_file3 = '/Users/jjb/Dropbox/ukidss_m_giants/data/mmt/final/GXAnd_f.fits'
sci_spec_file4 = '/Users/jjb/Dropbox/ukidss_m_giants/data/lbt/mstar3_spec.txt'
sci_spec_file5 = '/Users/jjb/Dropbox/ukidss_m_giants/data/lbt/mstar4_spec.txt'
sci_spec_file6 = '/Users/jjb/Dropbox/ukidss_m_giants/data/mmt/may2014/11May2014/final/Mg223.23_-0.81_f.fits'
sci_spec_file7 = '/Users/jjb/Dropbox/ukidss_m_giants/data/lbt/mstar10_spec.txt'
sci_spec_file8 = '/Users/jjb/Dropbox/ukidss_m_giants/data/lbt/mstar23_spec.txt'
sci_spec3 = read_mmt_spec(sci_spec_file3)
sci_spec4 = read_lbt_spec(sci_spec_file4)
sci_spec5 = read_lbt_spec(sci_spec_file5)
sci_spec6 = read_mmt_spec(sci_spec_file6)
sci_spec7 = read_lbt_spec(sci_spec_file7)
sci_spec8 = read_lbt_spec(sci_spec_file8)

smooth_flux = smooth(sci_spec.flux,smooth_beta)
smooth_flux2 = smooth(sci_spec2.flux,smooth_beta)
smooth_flux3 = smooth(sci_spec3.flux,smooth_beta)
smooth_flux4 = smooth(sci_spec4.flux,smooth_beta)
smooth_flux5 = smooth(sci_spec5.flux,smooth_beta)
smooth_flux6 = smooth(sci_spec6.flux,smooth_beta)
smooth_flux7 = smooth(sci_spec7.flux,smooth_beta)
smooth_flux8 = smooth(sci_spec8.flux,smooth_beta)

fig = plt.figure(figsize=(16, 10))
fig.subplots_adjust(bottom=0.05, top=0.95, left=0.15, right=0.95)
ax1 = fig.add_subplot(111)

ax1.plot(sci_spec.wave,smooth_flux)
ax1.plot(sci_spec2.wave,smooth_flux2+1)
ax1.plot(sci_spec6.wave,smooth_flux6+2)
ax1.plot(sci_spec3.wave,smooth_flux3+3)

ax1.plot(sci_spec4.wave,smooth_flux4+4)
ax1.plot(sci_spec5.wave,smooth_flux5+5)
ax1.plot(sci_spec7.wave,smooth_flux7+6)
#ax1.plot(sci_spec8.wave,smooth_flux8+7)

minorLocator   = MultipleLocator(20)
ax1.xaxis.set_minor_locator(minorLocator)
minorLocator   = MultipleLocator(0.2)
ax1.yaxis.set_minor_locator(minorLocator)

ax1.axvspan(8180, 8200, facecolor=(0.84313725,  0.89882353,  0.94431373))#, alpha=0.2)
ax1.axvspan(8493, 8501, facecolor=(0.97882353,  0.82039216,  0.82196078))
ax1.axvspan(8533, 8551, facecolor=(0.97882353,  0.82039216,  0.82196078))
ax1.axvspan(8652, 8669, facecolor=(0.97882353,  0.82039216,  0.82196078))
ax1.axvspan(8220, 8240, facecolor=(0.90039216,  0.90039216,  0.90039216))
ax1.axvspan(8293, 8300, facecolor=(0.90039216,  0.90039216,  0.90039216))
ax1.axvspan(8340, 8346, facecolor=(0.90039216,  0.90039216,  0.90039216))
ax1.axvspan(8394, 8402, facecolor=(0.90039216,  0.90039216,  0.90039216))
ax1.axvspan(8424, 8434, facecolor=(0.90039216,  0.90039216,  0.90039216))
ax1.axvspan(8459, 8466, facecolor=(0.90039216,  0.90039216,  0.90039216))



ax1.set_ylim(0.4,7.8)
ax1.set_xlim(8000,8800)
ax1.set_xlabel("Wavelength ($\AA$)", fontsize=35, style='italic')
ax1.set_ylabel("Normalized Flux + Constant", fontsize=35, style='italic')


#plt.show()
plt.savefig('/Users/jjb/Dropbox/ukidss_m_giants/figures/ukidss_lbt_spectra.png', dpi=300)
#plt.close()