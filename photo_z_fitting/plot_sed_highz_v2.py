"""
Routine to make MAGPHYS high-z (v2; Battisti et al. 2020) plots
Author: Andrew Battisti (Oct 2021), modified from code kindly provided by Katie Grasha and Tom Williams
"""
from __future__ import absolute_import, print_function, division
import os
import numpy as np
import astropy.units as u
from matplotlib.gridspec import GridSpec
#from astropy.cosmology import FlatLambdaCDM #MAGPHYS assumes FlatLambdaCDM, use same here for consistency
import matplotlib.pyplot as plt
#cosmo = FlatLambdaCDM(H0=70, Om0=0.3)


######  UPDATE THESE DIRECTORIES AND FILES  ######
# Directory where the MAGPHYS outputs are located and where figure is saved 
#os.chdir('/home/bijaya/Capstone_Research/MAGPHYS/magphys/magphys')

# Input user files
flux_file = '/Users/admin/data/observations.dat'
filter_file = '/Users/admin/data/filters.dat'

# Adjust as desired
sed_xrange=[1.01e-1,2e3]
sed_yrange=[8.01,13]
##################################################


# Read in galaxy list
galaxy_list = np.loadtxt(flux_file,dtype='str',usecols=0)
# Read in filter wavelengths
filter_wavelength = np.loadtxt(filter_file,usecols=1)

### Plot first few in the list###
#for i in range(3):
### Plot entire list###
for i in range(np.size(galaxy_list)):
    source = galaxy_list[i]
# Read in observed fluxes, errors, redshift
    obs_flux = np.loadtxt(source+'.fit',skiprows=2,max_rows=1)
    obs_err = np.loadtxt(source+'.fit',skiprows=3,max_rows=1)
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
    wavelength,attenuated_sed,unattenuated_sed = np.loadtxt(source+'.sed',skiprows=10,unpack=True)
    ar_inf = np.where(np.isinf(unattenuated_sed))
    unattenuated_sed[ar_inf] = np.nan

# Read in the PDFs
    Mstars_pdf = np.loadtxt(source+'.fit',skiprows=321,max_rows=70) #fine
    SFR_pdf = np.loadtxt(source+'.fit',skiprows=751,max_rows=60) #Fine
    sSFR_pdf = np.loadtxt(source+'.fit',skiprows=248,max_rows=70) #Fine
    Ldust_pdf = np.loadtxt(source+'.fit',skiprows=394,max_rows=70) #Fine
    Mdust_pdf = np.loadtxt(source+'.fit',skiprows=688,max_rows=60) #Fine
    age_M_pdf = np.loadtxt(source+'.fit',skiprows=814,max_rows=50) #Fine
    AV_pdf = np.loadtxt(source+'.fit',skiprows=1073,max_rows=80) #Fine
    #EbPrime_pdf = np.loadtxt(source+'.fit',skiprows=1223,max_rows=15) #Doesn't exist?
    mu_pdf = np.loadtxt(source+'.fit',skiprows=62,max_rows=20) #fine
    Tdust_pdf = np.loadtxt(source+'.fit',skiprows=1156,max_rows=14) #fine 

# Read in the percentiles
    Mstars_percentile = np.loadtxt(source+'.fit',skiprows=402,max_rows=1) #Fine
    SFR_percentile = np.loadtxt(source+'.fit',skiprows=812,max_rows=1) #Fine
    sSFR_percentile = np.loadtxt(source+'.fit',skiprows=319,max_rows=1) #Fine
    Ldust_percentile = np.loadtxt(source+'.fit',skiprows=465,max_rows=1) #Fine
    Mdust_percentile = np.loadtxt(source+'.fit',skiprows=749,max_rows=1) #Fine
    age_M_percentile = np.loadtxt(source+'.fit',skiprows=865,max_rows=1) #Fine
    AV_percentile = np.loadtxt(source+'.fit',skiprows=1154,max_rows=1) #Fine
    #EbPrime_percentile = np.loadtxt(source+'.fit',skiprows=1239,max_rows=1)
    mu_percentile = np.loadtxt(source+'.fit',skiprows=83,max_rows=1) #Fine
    Tdust_percentile = np.loadtxt(source+'.fit',skiprows=1171,max_rows=1) #fine

# Turn redshift to luminosity distance, and then fluxes to luminosities
#    Dlum = cosmo.luminosity_distance(redshift)
#    Dlum = Dlum.value

    x = 10**wavelength                                               
    y_at = np.log10(x*10**attenuated_sed) 
    y_un = np.log10(x*10**unattenuated_sed)
    
#convert to microns
    xx=x/1.e4

#==Lnu_sun->Llam_sun==
    L_flux=np.log10((1.+redshift)*obs_flux*3e14/filter_wavelength)

    L_eflux_lo=np.log10((1.+redshift)*obs_flux*3e14/filter_wavelength)- np.log10((1.+redshift)*obs_flux*3e14/filter_wavelength-(1.+redshift)*obs_err*3e14/filter_wavelength*1e4)    
    L_eflux_hi=-np.log10((1.+redshift)*obs_flux*3e14/filter_wavelength)+ np.log10((1.+redshift)*obs_flux*3e14/filter_wavelength+(1.+redshift)*obs_err*3e14/filter_wavelength)
  
    L_pflux=np.log10((1.+redshift)*magphys_flux*3e14/filter_wavelength)
    L_pflux_diff=np.log10((1.+redshift)*obs_flux*3e14/filter_wavelength)- np.log10((1.+redshift)*magphys_flux*3e14/filter_wavelength)

    if ar_limit: 
        L_flux_lim=np.log10((1.+redshift)*obs_flux_lim*3e14/filter_wavelength_lim)
        L_pflux_diff_lim=np.log10((1.+redshift)*obs_flux_lim*3e14/filter_wavelength_lim) -np.log10((1.+redshift)*magphys_flux_lim*3e14/filter_wavelength_lim)

    
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
    #ax10=fig.add_subplot(gs[3,2]) # Fourth row, third column
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



    ax3.plot(Mstars_pdf[:,0],Mstars_pdf[:,1], drawstyle='steps',c='black',lw=1)
    ax3.set_xlim([7.05,13.95])
    ax3.set_ylim([0,1])
    ax4.plot(SFR_pdf[:,0],SFR_pdf[:,1], drawstyle='steps',c='black',lw=1)
    ax4.set_xlim([-1.95,3.95])
    ax4.set_ylim([0,1])
    ax5.plot(sSFR_pdf[:,0],sSFR_pdf[:,1], drawstyle='steps',c='black',lw=1)
    ax5.set_xlim([-12.95,-6.05])
    ax5.set_ylim([0,1])
    ax6.plot(Ldust_pdf[:,0],Ldust_pdf[:,1], drawstyle='steps',c='black',lw=1)
    ax6.set_xlim([7.05,13.95])
    ax6.set_ylim([0,1])
    ax7.plot(Mdust_pdf[:,0],Mdust_pdf[:,1], drawstyle='steps',c='black',lw=1)
    ax7.set_xlim([4.05,9.95])
    ax7.set_ylim([0,1])
    ax8.plot(age_M_pdf[:,0],age_M_pdf[:,1], drawstyle='steps',c='black',lw=1)
    ax8.set_xlim([5.55,10.45])
    ax8.set_ylim([0,1])
    ax9.plot(AV_pdf[:,0],AV_pdf[:,1], drawstyle='steps',c='black',lw=1)
    ax9.set_xlim([0.01,19.875])
    ax9.set_ylim([0,1])
    #ax10.plot(EbPrime_pdf[:,0],EbPrime_pdf[:,1], drawstyle='steps',c='black',lw=1)
    #ax10.set_xlim([0.05,1.65])
    #ax10.set_ylim([0,1])
    ax11.plot(mu_pdf[:,0],mu_pdf[:,1], drawstyle='steps',c='black',lw=1)
    ax11.set_xlim([0.025,0.975])
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

    ax3.text(0.7, 0.875, '{:.2f}'.format(Mstars_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax3.transAxes, fontsize=12)
    ax3.text(0.97, 0.95, '+{:.2f}'.format(Mstars_percentile[3]-Mstars_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax3.transAxes, fontsize=12)
    ax3.text(0.97, 0.80, '-{:.2f}'.format(Mstars_percentile[2]-Mstars_percentile[1]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax3.transAxes, fontsize=12)
    ax4.text(0.7, 0.875, '{:.2f}'.format(SFR_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax4.transAxes, fontsize=12)
    ax4.text(0.97, 0.95, '+{:.2f}'.format(SFR_percentile[3]-SFR_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax4.transAxes, fontsize=12)
    ax4.text(0.97, 0.80, '-{:.2f}'.format(SFR_percentile[2]-SFR_percentile[1]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax4.transAxes, fontsize=12)
    ax5.text(0.7, 0.875, '{:.2f}'.format(sSFR_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax5.transAxes, fontsize=12)
    ax5.text(0.97, 0.95, '+{:.2f}'.format(sSFR_percentile[3]-sSFR_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax5.transAxes, fontsize=12)
    ax5.text(0.97, 0.80, '-{:.2f}'.format(sSFR_percentile[2]-sSFR_percentile[1]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax5.transAxes, fontsize=12)
    ax6.text(0.7, 0.875, '{:.2f}'.format(Ldust_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax6.transAxes, fontsize=12)
    ax6.text(0.97, 0.95, '+{:.2f}'.format(Ldust_percentile[3]-Ldust_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax6.transAxes, fontsize=12)
    ax6.text(0.97, 0.80, '-{:.2f}'.format(Ldust_percentile[2]-Ldust_percentile[1]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax6.transAxes, fontsize=12)
    ax7.text(0.7, 0.875, '{:.2f}'.format(Mdust_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax7.transAxes, fontsize=12)
    ax7.text(0.97, 0.95, '+{:.2f}'.format(Mdust_percentile[3]-Mdust_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax7.transAxes, fontsize=12)
    ax7.text(0.97, 0.80, '-{:.2f}'.format(Mdust_percentile[2]-Mdust_percentile[1]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax7.transAxes, fontsize=12)
    ax8.text(0.7, 0.875, '{:.2f}'.format(age_M_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax8.transAxes, fontsize=12)
    ax8.text(0.97, 0.95, '+{:.2f}'.format(age_M_percentile[3]-age_M_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax8.transAxes, fontsize=12)
    ax8.text(0.97, 0.80, '-{:.2f}'.format(age_M_percentile[2]-age_M_percentile[1]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax8.transAxes, fontsize=12)
    ax9.text(0.7, 0.875, '{:.2f}'.format(AV_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax9.transAxes, fontsize=12)
    ax9.text(0.97, 0.95, '+{:.2f}'.format(AV_percentile[3]-AV_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax9.transAxes, fontsize=12)
    ax9.text(0.97, 0.80, '-{:.2f}'.format(AV_percentile[2]-AV_percentile[1]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax9.transAxes, fontsize=12)
    #ax10.text(0.7, 0.875, '{:.2f}'.format(EbPrime_percentile[2]),
        #verticalalignment='top', horizontalalignment='right',
        #transform=ax10.transAxes, fontsize=12)
    #ax10.text(0.97, 0.95, '+{:.2f}'.format(EbPrime_percentile[3]-EbPrime_percentile[2]),
        #verticalalignment='top', horizontalalignment='right',
        #transform=ax10.transAxes, fontsize=12)
    #ax10.text(0.97, 0.80, '-{:.2f}'.format(EbPrime_percentile[2]-EbPrime_percentile[1]),
       #verticalalignment='top', horizontalalignment='right',
        #transform=ax10.transAxes, fontsize=12)
    ax11.text(0.7, 0.875, '{:.2f}'.format(mu_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax11.transAxes, fontsize=12)
    ax11.text(0.97, 0.95, '+{:.2f}'.format(mu_percentile[3]-mu_percentile[2]),
        verticalalignment='top', horizontalalignment='right',
        transform=ax11.transAxes, fontsize=12)
    ax11.text(0.97, 0.80, '-{:.2f}'.format(mu_percentile[2]-mu_percentile[1]),
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
    ax3.set_xlabel(r'$\log[M_*/M_\odot]$', fontsize=13)
    ax4.set_xlabel(r'$\log[$SFR$/(M_\odot$ yr$^{-1}$)]', fontsize=13)
    ax5.set_xlabel(r'$\log[$sSFR/yr$^{-1}$]', fontsize=13)
    ax6.set_xlabel(r'$\log[L_{dust}/L_\odot]$', fontsize=13)
    ax7.set_xlabel(r'$\log[M_{dust}/M_\odot]$', fontsize=13)

    ax8.set_ylabel('Likelihood', fontsize=14)
    ax8.set_xlabel(r'$\log[Age_M$/yr]', fontsize=13)
    ax9.set_xlabel(r'$A_V$/mag', fontsize=13)
    #ax10.set_xlabel(r"$E_b$'", fontsize=13)
    ax11.set_xlabel(r'$\mu$', fontsize=13)
    ax12.set_xlabel(r'$T_{dust}$/K', fontsize=13)

    plt.savefig(source+'.png')
    
print('Complete!')
