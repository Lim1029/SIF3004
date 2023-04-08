"""
Routine to make MAGPHYS+photo-z (Battisti et al. 2019) plots
Author: Andrew Battisti (Oct 2021), modified from code kindly provided by Katie Grasha and Tom Williams
"""
from __future__ import absolute_import, print_function, division
import os
import sys
import numpy as np
import astropy.units as u
from matplotlib.gridspec import GridSpec
from astropy.cosmology import FlatLambdaCDM #MAGPHYS assumes FlatLambdaCDM, use same here for consistency
import matplotlib.pyplot as plt
import pandas as pd
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)


######  UPDATE THESE DIRECTORIES AND FILES  ######

directory = sys.argv[1]
filter_path = sys.argv[2]
flux_path = sys.argv[3]

# Directory where the MAGPHYS outputs are located and where figure is saved 
os.chdir(directory)

# Input user files
flux_file = directory + flux_path 
filter_file = directory + filter_path 

# Adjust as desired
sed_xrange=[1.01e-1,3e5]
sed_yrange=[6.01,14]
##################################################


# Read in galaxy list
galaxy_list = np.loadtxt(flux_file,dtype='str',usecols=0)
galaxy_list = pd.read_table(flux_file,sep='\s+')['#name']
# print(galaxy_list)

# galaxy_list = galaxy_list[0:]
# print(galaxy_list)
# Read in filter wavelengths
filter_wavelength = np.loadtxt(filter_file,usecols=1)

### Plot first few in the list###
#for i in range(3):
### Plot entire list###
for i in range(len(galaxy_list)):
    source = galaxy_list[i]
    print(source)
# Read in observed fluxes, errors, photoz
    obs_flux = np.loadtxt(source+'.fit',skiprows=2,max_rows=1)
    obs_err = np.loadtxt(source+'.fit',skiprows=3,max_rows=1)
    print(obs_flux,obs_err)
    ar_limit = np.where(((obs_flux < 0) & (obs_err > 0)) | (obs_flux/obs_err < 1))
# Only do this if limits exist    
    if ar_limit: 
        filter_wavelength_lim=filter_wavelength[ar_limit]
        obs_flux_lim=obs_err[ar_limit]
   
    ar_neg = np.where((obs_flux < 0))
    obs_flux[ar_neg] = np.nan
    ar_neg = np.where((obs_err < 0))
    obs_err[ar_neg] = np.nan
    obs_flux[ar_limit] = np.nan
    obs_err[ar_limit] = np.nan
 
    tab = np.loadtxt(source+'.fit',skiprows=8,max_rows=1,unpack=True)
    chi2 = float(tab[2])
    redshift = float(tab[3])

    magphys_flux = np.loadtxt(source+'.fit',skiprows=12,max_rows=1)
    if ar_limit: 
        magphys_flux_lim=magphys_flux[ar_limit]
        
# Read in the SED
    wavelength,attenuated_sed,unattenuated_sed = np.loadtxt(source+'.sed',skiprows=11,unpack=True)
    ar_inf = np.where(np.isinf(unattenuated_sed))
    unattenuated_sed[ar_inf] = np.nan

# Read in the PDFs
    z_pdf = np.loadtxt(source+'.fit',skiprows=1144,max_rows=40)
    Mstars_pdf = np.loadtxt(source+'.fit',skiprows=241,max_rows=70)
    SFR_pdf = np.loadtxt(source+'.fit',skiprows=701,max_rows=60)
    sSFR_pdf = np.loadtxt(source+'.fit',skiprows=168,max_rows=70)
    Ldust_pdf = np.loadtxt(source+'.fit',skiprows=314,max_rows=70)
    Mdust_pdf = np.loadtxt(source+'.fit',skiprows=638,max_rows=60)
    age_M_pdf = np.loadtxt(source+'.fit',skiprows=764,max_rows=50)
    AV_pdf = np.loadtxt(source+'.fit',skiprows=1023,max_rows=80)
    EbPrime_pdf = np.loadtxt(source+'.fit',skiprows=1187,max_rows=15)
    Tdust_pdf = np.loadtxt(source+'.fit',skiprows=1106,max_rows=35)

# Read in the percentiles
    z_percentile = np.loadtxt(source+'.fit',skiprows=1185,max_rows=1)
    Mstars_percentile = np.loadtxt(source+'.fit',skiprows=312,max_rows=1)
    SFR_percentile = np.loadtxt(source+'.fit',skiprows=762,max_rows=1)
    sSFR_percentile = np.loadtxt(source+'.fit',skiprows=239,max_rows=1)
    Ldust_percentile = np.loadtxt(source+'.fit',skiprows=385,max_rows=1)
    Mdust_percentile = np.loadtxt(source+'.fit',skiprows=699,max_rows=1)
    age_M_percentile = np.loadtxt(source+'.fit',skiprows=815,max_rows=1)
    AV_percentile = np.loadtxt(source+'.fit',skiprows=1104,max_rows=1)
    EbPrime_percentile = np.loadtxt(source+'.fit',skiprows=1203,max_rows=1)
    Tdust_percentile = np.loadtxt(source+'.fit',skiprows=1142,max_rows=1)

# Turn redshift to luminosity distance, and then fluxes to luminosities
    Dlum = cosmo.luminosity_distance(redshift)
    Dlum = Dlum.value

#==1 Mpc = 3.08568e24 cm  ;4*!pi*3.08568e24^2/3.84e33 = 3.11588e16 ;Jy->Llam_sun==
    x = 10**wavelength                                               
    y_at = np.log10(10**attenuated_sed*3e-5/x*Dlum**2*3.11588e16) 
    y_un = np.log10(10**unattenuated_sed*3e-5/x*Dlum**2*3.11588e16)
    
#convert to microns
    xx=x/1.e4

    L_flux=np.log10(obs_flux*3e-5/(filter_wavelength*1e4)*Dlum**2*3.11588e16)

    L_eflux_lo=np.log10(obs_flux*3e-5/(filter_wavelength*1e4)*Dlum**2*3.11588e16)- np.log10(obs_flux*3e-5/(filter_wavelength*1e4)*Dlum**2*3.11588e16-obs_err*3e-5/(filter_wavelength*1e4)*Dlum**2*3.11588e16)
    L_eflux_hi=-np.log10(obs_flux*3e-5/(filter_wavelength*1e4)*Dlum**2*3.11588e16)+ np.log10(obs_flux*3e-5/(filter_wavelength*1e4)*Dlum**2*3.11588e16+obs_err*3e-5/(filter_wavelength*1e4)*Dlum**2*3.11588e16)

    L_pflux=np.log10(magphys_flux*3e-5/(filter_wavelength*1e4)*Dlum**2*3.11588e16)
    L_pflux_diff=np.log10(obs_flux*3e-5/(filter_wavelength*1e4)*Dlum**2*3.11588e16)- np.log10(magphys_flux*3e-5/(filter_wavelength*1e4)*Dlum**2*3.11588e16)

    if ar_limit: 
        L_flux_lim=np.log10(obs_flux_lim*3e-5/(filter_wavelength_lim*1e4)*Dlum**2*3.11588e16)
        L_pflux_diff_lim=np.log10(obs_flux_lim*3e-5/(filter_wavelength_lim*1e4)*Dlum**2*3.11588e16) -np.log10(magphys_flux_lim*3e-5/(filter_wavelength_lim*1e4)*Dlum**2*3.11588e16)

    
# Start plotting 
    fig = plt.figure(figsize=(10, 8)) #10 wide and 8 tall

    heights = [5, 1, 3, 3] # height of each row relative to each other

# have to set up two "gridspec", one for the top 2 rows that share x values
    gs1 = GridSpec(4,5, hspace=0, height_ratios=heights, bottom=0.2) # 4 rows, 5 columns

    ax1=fig.add_subplot(gs1[0,:]) # First row, span all columns
    ax2=fig.add_subplot(gs1[1,:], sharex=ax1) # Second row, span all columns


#second "gridspec" enviro, for the bottom 2 rows that share y values
    gs=GridSpec(4,5, wspace=0.0, height_ratios=heights) # 4 rows, 5 columns

    ax3=fig.add_subplot(gs[2,0]) # Third row, first column
    ax4=fig.add_subplot(gs[2,1]) # Third row, second column
    ax5=fig.add_subplot(gs[2,2]) # Third row, third column
    ax6=fig.add_subplot(gs[2,3]) # Third row, fourth column
    ax7=fig.add_subplot(gs[2,4]) # Third row, fifth column
    
    ax8=fig.add_subplot(gs[3,0]) # Fourth row, first column
    ax9=fig.add_subplot(gs[3,1]) # Fourth row, second column
    ax10=fig.add_subplot(gs[3,2]) # Fourth row, third column
    ax11=fig.add_subplot(gs[3,3]) # Fourth row, fourth column
    ax12=fig.add_subplot(gs[3,4]) # ThirdFourth row, fifth column


# remove the axis labels where x/y are shared
    plt.setp(ax1.get_xticklabels(), visible=False)

    plt.setp(ax4.get_yticklabels(), visible=False)
    plt.setp(ax5.get_yticklabels(), visible=False)
    plt.setp(ax6.get_yticklabels(), visible=False)
    plt.setp(ax7.get_yticklabels(), visible=False)
    
    plt.setp(ax9.get_yticklabels(), visible=False)
    plt.setp(ax10.get_yticklabels(), visible=False)
    plt.setp(ax11.get_yticklabels(), visible=False)
    plt.setp(ax12.get_yticklabels(), visible=False)
    


## plot your data
    ax1.plot(xx,y_un,c='cornflowerblue',lw=2)
    ax1.plot(xx,y_at,c='black',lw=2)

    ax1.plot(filter_wavelength,L_pflux,c='black',ls='none',marker='o', mfc='none',ms=10, mew=2)
    
    ax1.errorbar(filter_wavelength,L_flux,
             yerr=[L_eflux_lo,L_eflux_hi],c='red',
             ls='none',marker='o',zorder=20)
    if ar_limit: 
        ax1.errorbar(filter_wavelength_lim,L_flux_lim,yerr=0.5,uplims=1,c='red',zorder=20,linestyle='none')
        
    ax1.set_xlim(sed_xrange)
    ax1.set_ylim(sed_yrange)
    ax1.set_xscale('log')

    ax2.errorbar(filter_wavelength,L_pflux_diff,
                 yerr=[L_eflux_lo,L_eflux_hi],c='black',
                 ls='none',marker='o')
    if ar_limit: 
        ax2.errorbar(filter_wavelength_lim,L_pflux_diff_lim,yerr=0.5,uplims=1,c='black',zorder=20,linestyle='none')
#        for j in range(len(ar_limit)):
#            ax2.arrow(filter_wavelength_lim[j],L_pflux_diff_lim[j],0,1, head_width=3, head_length=6, fc='k', ec='k')
        
    ax2.axhline(y=0, color='black', linestyle='dashed',lw=1)
    ax2.set_ylim([-1,1])



    ax3.plot(z_pdf[:,0],z_pdf[:,1], drawstyle='steps',c='black',lw=1)
    ax3.set_xlim([0.0,7.9])
    ax3.set_ylim([0,1])
    ax4.plot(Mstars_pdf[:,0],Mstars_pdf[:,1], drawstyle='steps',c='black',lw=1)
    ax4.set_xlim([7.05,13.95])
    ax4.set_ylim([0,1])
    ax5.plot(SFR_pdf[:,0],SFR_pdf[:,1], drawstyle='steps',c='black',lw=1)
    ax5.set_xlim([-1.95,3.95])
    ax5.set_ylim([0,1])
    ax6.plot(sSFR_pdf[:,0],sSFR_pdf[:,1], drawstyle='steps',c='black',lw=1)
    ax6.set_xlim([-12.95,-6.05])
    ax6.set_ylim([0,1])
    ax7.plot(Ldust_pdf[:,0],Ldust_pdf[:,1], drawstyle='steps',c='black',lw=1)
    ax7.set_xlim([7.05,13.95])
    ax7.set_ylim([0,1])
    ax8.plot(Mdust_pdf[:,0],Mdust_pdf[:,1], drawstyle='steps',c='black',lw=1)
    ax8.set_xlim([6.05,11.95])
    ax8.set_ylim([0,1])
    ax9.plot(age_M_pdf[:,0],age_M_pdf[:,1], drawstyle='steps',c='black',lw=1)
    ax9.set_xlim([5.55,10.45])
    ax9.set_ylim([0,1])
    ax10.plot(AV_pdf[:,0],AV_pdf[:,1], drawstyle='steps',c='black',lw=1)
    ax10.set_xlim([0.01,19.875])
    ax10.set_ylim([0,1])
    ax11.plot(EbPrime_pdf[:,0],EbPrime_pdf[:,1], drawstyle='steps',c='black',lw=1)
    ax11.set_xlim([0.05,1.65])
    ax11.set_ylim([0,1])
    ax12.plot(Tdust_pdf[:,0],Tdust_pdf[:,1], drawstyle='steps',c='black',lw=1)
    ax12.set_xlim([11,80])
    ax12.set_ylim([0,1])
    
# Output results
    ax1.text(0.02, 0.95, 'ID: '+source,
        verticalalignment='top', horizontalalignment='left',
        transform=ax1.transAxes, fontsize=15)
    ax1.text(0.98, 0.95, r'$z_{fit}=$'+'{:.2f}'.format(redshift),
        verticalalignment='top', horizontalalignment='right',
        transform=ax1.transAxes, fontsize=15)
    ax1.text(0.98, 0.85, r'$\chi^2=$'+'{:.2f}'.format(chi2),
        verticalalignment='top', horizontalalignment='right',
        transform=ax1.transAxes, fontsize=15)

    ax3.text(0.7, 0.875, '{:.2f}'.format(z_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax3.transAxes, fontsize=12)
    ax3.text(0.97, 0.95, '+{:.2f}'.format(z_percentile[3]-z_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax3.transAxes, fontsize=12)
    ax3.text(0.97, 0.80, '-{:.2f}'.format(z_percentile[2]-z_percentile[1]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax3.transAxes, fontsize=12)
    ax4.text(0.7, 0.875, '{:.2f}'.format(Mstars_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax4.transAxes, fontsize=12)
    ax4.text(0.97, 0.95, '+{:.2f}'.format(Mstars_percentile[3]-Mstars_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax4.transAxes, fontsize=12)
    ax4.text(0.97, 0.80, '-{:.2f}'.format(Mstars_percentile[2]-Mstars_percentile[1]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax4.transAxes, fontsize=12)
    ax5.text(0.7, 0.875, '{:.2f}'.format(SFR_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax5.transAxes, fontsize=12)
    ax5.text(0.97, 0.95, '+{:.2f}'.format(SFR_percentile[3]-SFR_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax5.transAxes, fontsize=12)
    ax5.text(0.97, 0.80, '-{:.2f}'.format(SFR_percentile[2]-SFR_percentile[1]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax5.transAxes, fontsize=12)
    ax6.text(0.7, 0.875, '{:.2f}'.format(sSFR_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax6.transAxes, fontsize=12)
    ax6.text(0.97, 0.95, '+{:.2f}'.format(sSFR_percentile[3]-sSFR_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax6.transAxes, fontsize=12)
    ax6.text(0.97, 0.80, '-{:.2f}'.format(sSFR_percentile[2]-sSFR_percentile[1]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax6.transAxes, fontsize=12)
    ax7.text(0.7, 0.875, '{:.2f}'.format(Ldust_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax7.transAxes, fontsize=12)
    ax7.text(0.97, 0.95, '+{:.2f}'.format(Ldust_percentile[3]-Ldust_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax7.transAxes, fontsize=12)
    ax7.text(0.97, 0.80, '-{:.2f}'.format(Ldust_percentile[2]-Ldust_percentile[1]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax7.transAxes, fontsize=12)
    ax8.text(0.7, 0.875, '{:.2f}'.format(Mdust_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax8.transAxes, fontsize=12)
    ax8.text(0.97, 0.95, '+{:.2f}'.format(Mdust_percentile[3]-Mdust_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax8.transAxes, fontsize=12)
    ax8.text(0.97, 0.80, '-{:.2f}'.format(Mdust_percentile[2]-Mdust_percentile[1]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax8.transAxes, fontsize=12)
    ax9.text(0.7, 0.875, '{:.2f}'.format(age_M_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax9.transAxes, fontsize=12)
    ax9.text(0.97, 0.95, '+{:.2f}'.format(age_M_percentile[3]-age_M_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax9.transAxes, fontsize=12)
    ax9.text(0.97, 0.80, '-{:.2f}'.format(age_M_percentile[2]-age_M_percentile[1]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax9.transAxes, fontsize=12)
    ax10.text(0.7, 0.875, '{:.2f}'.format(AV_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax10.transAxes, fontsize=12)
    ax10.text(0.97, 0.95, '+{:.2f}'.format(AV_percentile[3]-AV_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax10.transAxes, fontsize=12)
    ax10.text(0.97, 0.80, '-{:.2f}'.format(AV_percentile[2]-AV_percentile[1]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax10.transAxes, fontsize=12)
    ax11.text(0.7, 0.875, '{:.2f}'.format(EbPrime_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax11.transAxes, fontsize=12)
    ax11.text(0.97, 0.95, '+{:.2f}'.format(EbPrime_percentile[3]-EbPrime_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax11.transAxes, fontsize=12)
    ax11.text(0.97, 0.80, '-{:.2f}'.format(EbPrime_percentile[2]-EbPrime_percentile[1]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax11.transAxes, fontsize=12)
    ax12.text(0.7, 0.875, '{:.2f}'.format(Tdust_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax12.transAxes, fontsize=12)
    ax12.text(0.97, 0.95, '+{:.2f}'.format(Tdust_percentile[3]-Tdust_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax12.transAxes, fontsize=12)
    ax12.text(0.97, 0.80, '-{:.2f}'.format(Tdust_percentile[2]-Tdust_percentile[1]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax12.transAxes, fontsize=12)

# set the x/y labels
    ax1.set_ylabel(r'$\log(\lambda L_\lambda$/$L_\odot$)', fontsize=18)
    ax2.set_ylabel('Resid.', fontsize=18)
    ax2.set_xlabel(r'$\lambda$/$\mu$m [obs-frame]', fontsize=18)

    ax3.set_ylabel('Likelihood', fontsize=14)
    ax3.set_xlabel(r'$z_{phot}$', fontsize=13)
    ax4.set_xlabel(r'$\log[M_*/M_\odot]$', fontsize=13)
    ax5.set_xlabel(r'$\log[$SFR$/(M_\odot$ yr$^{-1}$)]', fontsize=13)
    ax6.set_xlabel(r'$\log[$sSFR/yr$^{-1}$]', fontsize=13)
    ax7.set_xlabel(r'$\log[L_{dust}/L_\odot]$', fontsize=13)

    ax8.set_ylabel('Likelihood', fontsize=14)
    ax8.set_xlabel(r'$\log[M_{dust}/M_\odot]$', fontsize=13)
    ax9.set_xlabel(r'$\log[Age_M$/yr]', fontsize=13)
    ax10.set_xlabel(r'$A_V$/mag', fontsize=13)
    ax11.set_xlabel(r"$E_b$'", fontsize=13)
    ax12.set_xlabel(r'$T_{dust}$/K', fontsize=13)

    plt.savefig(source+'.png')
    
print('Complete!')
