#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Charles Fortenbach

Code associated with paper submitted to AAS Journals (6/10/19)

Title: A FRAMEWORK FOR OPTIMIZING EXOPLANET TARGET SELECTION FOR
THE JAMES WEBB SPACE TELESCOPE

Subprogram:  ExoT_Master.py

Version:  1.0

"""

import numpy as np
#import mr_forecast as mr
import subprocess
import gc
import ephem
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.interpolate import UnivariateSpline
import os



###############################################################################
# Read in JET input file
###############################################################################
JET_params = np.genfromtxt('JET_Input.txt', dtype=str, delimiter=', ')

#JWST instrument/mode
Instrument = JET_params[0, 1]

# Instrument noise floor (ppm)
nfloor = float(JET_params[1, 1])

# Set common plot parameters (e.g., for NIRSpec G235M):
wavelimlo = float(JET_params[2, 1])   #  low end of plot wavelength scale (nm)
wavelimhi = float(JET_params[3, 1])   # high end of plot wavelength scale (nm)

# Set equation of state for model atmospheres, low and high metallicity
EOS_lo = JET_params[4, 1]
EOS_hi = JET_params[5, 1]

# pressure of cloud top (in Pa), -- or leave at 0.0 for no clouds
Cloud_lo = float(JET_params[6, 1])
Cloud_hi = float(JET_params[7, 1])

# Spectral resolution
Res = int(JET_params[8, 1])

# Saturation magnitude
Jmag_sat = float(JET_params[9, 1])

# Set starting and ending targets from survey list
Nrow_start = int(JET_params[10, 1])
Nrow_end = int(JET_params[11, 1])

# Set # of dBIC sample runs for each transit case
Npdxo = int(JET_params[12, 1])

# Set dBIC threshold
dBIC_thresh = int(JET_params[13, 1])

# Set Run flags
RunExoT = JET_params[14, 1]
RunPdxo = JET_params[15, 1]
RunRank = JET_params[16, 1]



###############################################################################
# Check Run flag
###############################################################################
if RunExoT == 'N':
    raise SystemExit
    
    
    
###############################################################################
# Check for target list in sequence
###############################################################################
summary_data_old = np.loadtxt('Pdxo_Output.txt', ndmin=2, skiprows=0)
targets_old = np.max(summary_data_old, axis=0)
last_target_old = targets_old[2]

if Nrow_start != last_target_old + 1:
    print('\n')
    print(' !Warning: The starting target value (Nrow_start) is '
          'not in sequence.')
    print('\n')
    print(' If you proceed you may overwrite the previous'
          ' entries.')
    print('\n')
    xx = input(' Proceed with ExoT? (Y/N): ')
    print('\n')
    if xx != 'Y':
        print ('   exiting ExoT . . .')
        raise SystemExit
        


###############################################################################
# Print notification of running ExoT_Master
###############################################################################
print('\n')                      
print(" Running ExoT_Master:")



###############################################################################
# Read in survey data file 
###############################################################################
print('\n') 
print(" Reading in survey data . . . " ) 
survey_data = np.loadtxt('target_survey.txt', skiprows=32)



###############################################################################
# Initialize data transfer array
###############################################################################
Nrow = Nrow_end - Nrow_start + 1
Ncol = 15
transfer_data = np.zeros((Nrow, Ncol))



###############################################################################
# Set constants
###############################################################################
Nsample_size = 12500    # Chen & Kipping R-M est. parameters (MCMC)
Ngrid_size = 12500      # Chen & Kipping R-M est. parameters (MCMC)
Rfactr = 0.12           # Rp_stdev / Rp_mean guestimate

G = 6.6742867e-11       # Gravitational const., N m^2/kg^2
Mearth = 5.9736e24      # Earth mass, kg
Rearth = 6.378136e06    # Earth radius, m
AU = 1.4959787066e11    # AU, m
Rsun = 6.95508e08       # Solar radius, m
Teffsun = 5777          # effective temp of Sun, K



###############################################################################
#  Compute model atmospheric transmission spectra   
###############################################################################
print('\n') 
print(" Computing model transmission spectra using Exo-Transmit:\n" ) 


# Loop on target index 
for i in range(Nrow_start-1, Nrow_end):

    # Set display counter, n (starts at 1, while i starts at 0)
    n = i + 1
    
    #  Print progress note
    print("      Generating model spectra for target: ", str(n), "\n") 

    ###########################################################################  
    # Extract survey data for a given target
    ###########################################################################
    RAdeg = survey_data[i, 0]     # decimal RA (equatorial)
    DEdeg = survey_data[i, 1]     # decimal Dec (equatorial)
    Rp_mean = survey_data[i, 2]   # mean radius of planet in Earth radii
    Rp_std = Rp_mean * Rfactr     # stdev of radius of planet in Earth radii
    Per = survey_data[i, 3]       # orbital period in days
    S = survey_data[i, 4]         # insolation in Earth units
    K = survey_data[i, 5]         # radial vel. semi-amplitude (m/s)
    Rstar = survey_data[i, 6]     # stellar radius in Solar radii
    Teff = survey_data[i, 7]      # target star's effective temp in K
    Jmag = survey_data[i, 10]     # Jmag of target star
    
    
    ###########################################################################  
    # Compute additional planetary parameters including planet mass   
    ########################################################################### 
    
    # Compute semi-major axis of planet's orbit (assuming circular orbits)
    aAU = (1/np.sqrt(S))*(Rstar)*(Teff/Teffsun)**2     # in AU
    a = aAU*(AU/Rsun)                                  # in Solar radii
    
    # Compute planet's equilibrium temp. (K)
    Teq = Teff*np.sqrt(Rstar/(2*a))           
   
    
    # Estimate planet mass
    
    ############################# Inactive ####################################
    #
    # Compute planet mass using Chen Kipping method  (earth masses)
    #
    # Mmed, Mplus, Mminus = mr.Rstat2M(mean=Rp_mean, std=Rp_std, unit='Earth', 
    #                                 sample_size=Nsample_size, 
    #                                 grid_size=Ngrid_size) 
    #
    ###########################################################################
    
    # Forecaster power laws (see Chen & Kipping 2016)
    if Rp_mean <= 1.23:
        cc1 = 2.04/(1.23**(1/0.279))            
        Mmed = cc1*Rp_mean**(1/0.279)
    
    if Rp_mean > 1.23:
        cc2 = 2.04/(1.23**(1/0.589))            
        Mmed = cc2*Rp_mean**(1/0.589)
        
    Mplus = 0.     # only interested in mean value of planet mass
    Mminus = 0.    # only interested in mean value of planet mass

    # Compute planet's surface gravity (m/s^2)
    gs = G*Mmed*Mearth/((Rp_mean*Rearth)**2)  
    
    # Compute median planet density (earth units)
    Rhomed = Mmed/(Rp_mean**3)                
    
    # Compute transit duration in hrs (circular orbits, and i = 90 deg)
    # tdur = (1 + (Rp_mean*Rearth/Rsun)/Rstar)*(13)*(Per/365)*(Rstar/aAU)
    tdur = (1 + (Rp_mean*Rearth/Rsun)/Rstar)*(Rstar*Per*24)/(np.pi*a)


    # JWST mission days:  10yr * 365 d/yr less 6mo. commissioning
    tmission = (10-0.5)*365
               
    # transits available (not necessarily observable) in 10 yr mission           
    nt10yr_available = (1/Per)*tmission       
    
    # Convert decimal coords to hms
    c = SkyCoord(ra=RAdeg*u.degree, dec=DEdeg*u.degree, frame='icrs')
    
    cdms_ra  = (str(c.ra.hms[0]) + ':' + str(c.ra.hms[1]) + ':' 
                + str(c.ra.hms[2]))
    cdms_dec = (str(c.dec.dms[0]) + ':'+str(np.abs(c.dec.dms[1])) + ':'
                + str(np.abs(c.dec.dms[2])))
    
    # Convert to Ecliptic coordinates 
    m = ephem.Equatorial(cdms_ra, cdms_dec)
    ecl = ephem.Ecliptic(m)
    
    # Change variable type
    zzlat = float(ecl.lat)
    
    # Convert to degrees and round
    zzlat_deg = round(zzlat*180/np.pi, 1)
    
    # For this ecliptic latitude, estimate observable days per year
    if np.abs(zzlat_deg) <= 45:
        eclat_lo = [0, 5, 10, 20, 30, 35, 40, 43, 45]
        obsdays_lo = [100, 101, 103, 110, 125, 136, 154, 175, 200]
        z_lo = UnivariateSpline(eclat_lo, obsdays_lo, s=1)
        obsdays = np.round(z_lo(np.abs(zzlat_deg)),0)
    
    elif np.abs(zzlat_deg) > 45 and np.abs(zzlat_deg) < 85:
        eclat_hi = [45, 50, 60, 70, 75, 80, 85]
        obsdays_hi = [200, 201, 204, 210, 225, 255, 360]
        z_hi = UnivariateSpline(eclat_hi, obsdays_hi, s=1)
        obsdays = np.round(z_hi(np.abs(zzlat_deg)),0)
        
    elif np.abs(zzlat_deg) >= 85:
        obsdays = 365
    
    # For how much of year is target observable
    obsdays_factr = np.round(obsdays/365, 2)
      
    # Compute number of transits observable in 10 yr mission
    # Assumes no overlapping transits
    nt10yr = nt10yr_available*obsdays_factr              
    
    
    ###########################################################################  
    # Determine target planet's demographic category
    ###########################################################################  
    if Teq < 400 and Rp_mean < 1.7:
        Cat = 1
                 
    elif Teq >= 400 and Teq <= 800 and Rp_mean < 1.7: 
        Cat = 2
        
    elif Teq > 800 and Rp_mean < 1.7:    
        Cat = 3
               
    elif Teq < 400 and Rp_mean >= 1.7 and Rp_mean <= 4.0:
        Cat = 4
             
    elif Teq >= 400 and Teq <= 800 and Rp_mean >= 1.7 and Rp_mean <= 4.0:
        Cat = 5
             
    elif Teq > 800 and Rp_mean >= 1.7 and Rp_mean <= 4.0:
        Cat = 6
        
    elif Rp_mean > 4.0:
        Cat = 7        
    
    
    ###########################################################################  
    # Define transfer data for each target
    ###########################################################################    
    ntransfer = i-(Nrow_start-1)
    
    transfer_data[ntransfer, 0]  = Nrow_start + ntransfer
    transfer_data[ntransfer, 1]  = a
    transfer_data[ntransfer, 2]  = Teq
    transfer_data[ntransfer, 3]  = Mmed
    transfer_data[ntransfer, 4]  = Mplus
    transfer_data[ntransfer, 5]  = Mminus
    transfer_data[ntransfer, 6]  = gs
    transfer_data[ntransfer, 7]  = Rhomed
    transfer_data[ntransfer, 8]  = tdur
    transfer_data[ntransfer, 9]  = nt10yr    
    transfer_data[ntransfer, 10] = Cat


    # Loop on two generic atmosphere models (low and high metallicity)    
    for j in range(0, 2):        
    
    
        #######################################################################
        # Set up userInput.in file for Exo_Transmit
        #######################################################################
            
        fout = open("userInput.in", "w", newline='')
        fout.write("userInput.in - \n")
        fout.write(("Formatting here is very important, please put the "
                    "instructed content on the instructed line.\n"))
        fout.write("Exo_Transmit home directory:\n")
        fout.write(os.getcwd() + "\n")
            
        # Determine target planet's P-T profile based on Teq  
        t_P = int(round(Teq/100, 0)*100)
        
        if t_P < 300:
            t_P = 300
            
        elif t_P > 1500:
            t_P = 1500
            
        t_P_filename = "/T_P/t_p_" + str(t_P) + "K.dat"
        fout.write("Temperature-Pressure data file:\n")
        fout.write(t_P_filename)
            
        # Select the atmospheric equation of state file (EOS)    
        if j == 0:        
            EOS = EOS_lo
        
        elif j == 1:
            EOS = EOS_hi 
            
        EOS_filename =  "/EOS/" + str(EOS) + ".dat"
        fout.write("\nEquation of State file:\n")
        fout.write(EOS_filename)
            
        # Enter Output Spectrum filename (e.g. Test_Spec1.txt) 
        if j == 0:        
            Spectrum_filename = ("/Spectra/" + "Trans_Spec_ExoT_" 
                                 + str(n) + "_lo_metal.txt")
              
        elif j == 1:
            Spectrum_filename = ("/Spectra/" + "Trans_Spec_ExoT_" 
                                 + str(n) + "_hi_metal.txt") 
        
        fout.write("\nOutput file:\n")
        fout.write(Spectrum_filename)
        
        # Enter planet surface gravity, gs (in m/s^-2)
        gs_in = str(gs)
        fout.write("\nPlanet surface gravity (in m/s^-2):\n")   
        fout.write(gs_in)
              
        # Enter planet radius, Rp (in m): 
        # -- This is the planet's radius at the base of the atmosphere 
        # (or at the cloud top for cloudy calculations
        Rp_in = str(Rp_mean*Rearth)
        fout.write(("\nPlanet radius (in m): -- This is the planet's radius at"
                    " the base of the atmosphere (or at the cloud top for"
                    " cloudy calculations)\n"))
        fout.write(Rp_in)
                  
        # Enter star radius, Rstar (m)
        Rstar_in = str(Rstar*Rsun)
        fout.write("\nStar radius (in m):\n")
        fout.write(Rstar_in)
                      
        # Select the pressure of cloud top (in Pa), 
        # -- or leave at 0.0 if you want no cloud calculations   
        if j == 0:        
            Pcloud_in = Cloud_lo
        
        elif j == 1:
            Pcloud_in = Cloud_hi 
            
        Pcloud_in = str(Pcloud_in) + ".0"
        fout.write(("\nPressure of cloud top (in Pa), "
                    "-- or leave at 0.0 if you want no cloud calculations.\n"))
        fout.write(Pcloud_in)
                       
        # Enter Rayleigh scattering augmentation factor -- default is 1.0.  
        # Can be increased to simulate additional sources of scattering.  
        # 0.0 turns off Rayleigh scattering.
        Rayleigh_in = '1'
        Rayleigh_in = str(Rayleigh_in) + ".0" 
        fout.write(("\nRayleigh scattering augmentation factor"
                    " -- default is 1.0.  Can be increased to simulate "
                    "additional sources of scattering.  "
                    "0.0 turns off Rayleigh scattering.\n"))
        fout.write(Rayleigh_in)
                         
        fout.write("\nEnd of userInput.in (Do not change this line)")
        fout.close()
        
           
        #######################################################################
        # Set up selectChem.in file for Exo_Transmit - manual input
        #######################################################################
        #print(' ')
        #print('Enter data to build Exo_Transmit selectChem.in file: \n')
        #print('     ==>> Currently using default values')
            
                          
        #######################################################################
        # Call Exo_Transmit executable and generate transmission spectrum
        #######################################################################  
        gc.collect()
        p = subprocess.Popen(os.getcwd() + '/./Exo_Transmit', 
                             stdout=subprocess.PIPE, shell=True)
        (output, err) = p.communicate()
            


###############################################################################
# Write transfer_data to Pdxo input file
###############################################################################
np.savetxt('Pdxo_Input.txt', transfer_data)

raise SystemExit