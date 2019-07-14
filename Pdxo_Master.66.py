

"""
@author: Charles Fortenbach

Code associated with paper submitted to AAS Journals (6/10/19)

Title: A FRAMEWORK FOR OPTIMIZING EXOPLANET TARGET SELECTION FOR
THE JAMES WEBB SPACE TELESCOPE

Subprogram:  Pdxo_Master.py

Version:  1.0

This project is licensed under the GNU GPLv3 License.

"""

import warnings
warnings.filterwarnings('ignore')
import pandexo.engine.justdoit as jdi 
import pandexo.engine.justplotit as jpi
import pickle as pk
import numpy as np
import matplotlib.pyplot as plt
import os
from spectres import spectres
from scipy import optimize
from scipy.interpolate import UnivariateSpline



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
if RunPdxo == 'N':
    raise SystemExit
    
    
    
###############################################################################
# Check sequential targets
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
    xx = input(' Proceed with Pdxo? (Y/N): ')
    print('\n')
    if xx != 'Y':
        print ('   exiting Pdxo . . .')
        raise SystemExit
    
    
###############################################################################
# Print notification of running Pdxo_Master
###############################################################################
print('\n')                      
print(" Running Pdxo_Master:")



###############################################################################
# Read in survey data file 
###############################################################################
survey_data = np.loadtxt('target_survey.txt', skiprows=32)



###############################################################################
# Initialize transit cases
###############################################################################
ints = [1, 2, 4, 6, 8, 10, 13, 16, 20, 25, 30, 35, 40, 45, 50]



###############################################################################
# Initialize arrays
###############################################################################
Nrow = Nrow_end - Nrow_start + 1
nt_lo = np.zeros(Nrow, dtype=int)
nt_hi = np.zeros(Nrow, dtype=int)
tT_lo = np.zeros(Nrow, dtype=float)
tT_hi = np.zeros(Nrow, dtype=float)
tT_avg = np.zeros(Nrow, dtype=float)



###############################################################################
# Transfer planetary system data to Pdxo_Master 
###############################################################################
print('\n') 
print(" Transferring planetary system data and model spectra to Pandexo . . .")
pdxinput = np.loadtxt('Pdxo_Input.txt', ndmin=2, skiprows=0)



###############################################################################
#  Compute transits needed for detection given instrument noise characteristics
###############################################################################

# Loop on target index
for i in range(Nrow_start-1, Nrow_end):

    # Set display counter, n (starts at 1, while i starts at 0)
    n = i + 1
    
    # Set pdxinput counter          
    ntransfer = n - int(pdxinput[0, 0])
    
    # Print progress note
    print('\n\n', '     Generating simulated instrument spectra '
           'for target: ', str(n),'\n')
    
    
    
    # We are only considering planet radius < 10 Rearth    
    Rp = survey_data[i, 2]
    if Rp < 10:
      
    
        # Instrument has mag threshold:  
        # (e.g., nothing brighter than ~ Jmag 7.8 at 2500K  for NIRSpec G235M)  
        Jmag = survey_data[i, 10]
        if Jmag > Jmag_sat:
    
            
            # Loop on two generic atmosphere models (low and high metallicity)    
            for j in range(0, 2): 
                
                ###############################################################
                # Reset arrays
                ###############################################################
                xnt = np.zeros((len(ints), Npdxo), dtype=int)
                dBIC = np.zeros((len(ints), Npdxo), dtype=float)
                dBIC_mean = np.zeros((len(ints)), dtype=float)
                dBIC_sdev = np.zeros((len(ints)), dtype=float)
                xx = np.zeros(len(ints), dtype=float)
                yy = np.zeros(len(ints), dtype=float)
                
        
                # Read model transmission spectrum data file from working dir.
                if j == 0:
                    std_spectra = np.loadtxt(os.getcwd() + '/Spectra/'
                                             + 'Trans_Spec_ExoT_'+ str(n) 
                                             + '_lo_metal.txt', skiprows = 2)
                    
                elif j == 1:
                    std_spectra = np.loadtxt(os.getcwd() + '/Spectra/'
                                             + 'Trans_Spec_ExoT_'+ str(n) 
                                             + '_hi_metal.txt', skiprows = 2)
    
                                          
                # Match model spectrum limits to instrument bandpass
                Ndata_low = np.argmin(np.abs(std_spectra[:, 0]*10**6 
                                             - wavelimlo)) - 4
                Ndata_hi = np.argmin(np.abs(std_spectra[:, 0]*10**6 
                                            - wavelimhi)) + 8
                Nspread = Ndata_hi - Ndata_low
                std_spectra_1 = np.zeros((Nspread, 2), dtype=float)
                
                for jj in range(0, Nspread):
                    # 100 cm per meter
                    std_spectra_1[jj, 0] = std_spectra[jj + Ndata_low, 0]*100 
                    # convert % to actual
                    std_spectra_1[jj, 1] = std_spectra[jj + Ndata_low, 1]/100   
    
    
                # Save matched model spectrum to working directory    
                if j == 0:    
                    np.savetxt(os.getcwd() + '/Spectra/' 
                               + 'Trans_Spec_Xfer_' + str(n) + '_lo_metal.txt', 
                               std_spectra_1, fmt="%12.10f  %12.10f")
                
                elif j == 1:
                    np.savetxt(os.getcwd() + '/Spectra/' 
                               + 'Trans_Spec_Xfer_' + str(n) + '_hi_metal.txt', 
                               std_spectra_1, fmt="%12.10f  %12.10f")
    
                            
                # Set up input dictionaries for Pandexo
                exo_dict = jdi.load_exo_dict()
                
                # saturation level in percent of full well
                exo_dict['observation']['sat_level'] = 80              
                exo_dict['observation']['sat_unit'] = '%'
                
                # fixed binning. 
                exo_dict['observation']['R'] = Res
                
                # defines how you specify out of transit observing time
                # 'frac' : fraction of time in transit versus out = in/out                                                           
                exo_dict['observation']['baseline_unit'] = 'total'    
                                                                                                                                       
                total_obs = 2*pdxinput[ntransfer, 8]*3600
                
                # in accordance with what was specified above 
                # (total observing time in sec)                 
                exo_dict['observation']['baseline'] = total_obs
                
                # this can be a fixed level or it can be a filepath (ppm)
                exo_dict['observation']['noise_floor'] = nfloor       
            
                # phoenix or user (if you have your own) 
                exo_dict['star']['type'] = 'phoenix'

                # magnitude of the system                      
                exo_dict['star']['mag'] = Jmag   
                
                # for J mag = 1.25, H = 1.6, K = 2.22.. etc (all in microns)                       
                exo_dict['star']['ref_wave'] = 1.25        
                
                # Teff in K 
                exo_dict['star']['temp'] = survey_data[i, 7] 
                
                # as log Fe/H
                exo_dict['star']['metal'] = 0.0                       
                exo_dict['star']['logg'] = 4.0
                
                # tells pandexo you are uploading your own spectrum
                exo_dict['planet']['type'] ='user' 

                # options include "um","nm" ,"Angs", "sec" (for phase curves)                   
                exo_dict['planet']['w_unit'] = 'cm'    
                
                # other options are 'fp/f*' , 'rp^2/r*^2'
                exo_dict['planet']['f_unit'] = 'rp^2/r*^2'            
                
                tdur_sec = pdxinput[ntransfer, 8]*3600 
                
                # transit duration                                                
                exo_dict['planet']['transit_duration'] = tdur_sec
                
                # any unit of time in accordance with astropy.units
                exo_dict['planet']['td_unit'] = 's'                   
    
    
                # Define model spectrum for Pandexo dictionary
                if j == 0:    
                    exo_dict['planet']['exopath'] = (os.getcwd() + '/Spectra/' 
                            + 'Trans_Spec_Xfer_' + str(n) + '_lo_metal.txt')
                    print('\n', '        Comparing simulation to model spectrum'
                           ' and flat line for lo_metal atm:','\n')
                    
                
                elif j == 1:
                    exo_dict['planet']['exopath'] = (os.getcwd() + '/Spectra/' 
                            + 'Trans_Spec_Xfer_' + str(n) + '_hi_metal.txt')
                    print('\n', '        Comparing simulation to model spectrum'
                           ' and flat line for hi_metal atm:','\n')
                 
                   
                ###############################################################
                # Run Pandexo for a single transit
                ###############################################################             
    
                # Set for single transit             
                exo_dict['observation']['noccultations'] = 1   
                                
                # Run Pandexo
                single_transit = jdi.run_pandexo(exo_dict, [Instrument], 
                                                 save_file=True, 
                                                 output_path=os.getcwd(), 
                                                 output_file = 'singlerun.p')
                       
                # Extract output from pickle file
                out = pk.load(open('singlerun.p','rb'))
                
                # Generate simulated instrument spectrum with errors
                x, y, e = jpi.jwst_1d_spec(out, model=True, 
                                           x_range=[wavelimlo, wavelimhi], 
                                           plot=False)
    
                
                ###############################################################
                # Plot simulated instr. spectrum with model spectrum overlay
                ############################################################### 
                x1 = np.ndarray.flatten(np.array(x))
                y100 = np.ndarray.flatten(np.asarray(y)*100)
                e100 = np.ndarray.flatten(np.asarray(e)*100)
                
                                   
                # Rebin xfer data from Exo-Transmit to Pandexo output res
                # Call the spectres function to resample the input spectrum 
                # or spectra to the new wavelength grid
                rebinned_model = spectres(x1, std_spectra_1[:, 0]*10**4, 
                                          std_spectra_1[:, 1])               
                y_rebin = rebinned_model*100

                fignum = 1
                fname1 = [os.getcwd() + '/Spectra/' 
                          + 'Trans_Spec_Pdxo_' + str(n) + '_lo_metal.svg', 
                          os.getcwd() + '/Spectra/' 
                          + 'Trans_Spec_Pdxo_' + str(n) + '_hi_metal.svg']           
                
                plt.figure(fignum, figsize = (17, 9))
                plt.xlabel(r'Wavelength ($\mu$m)', fontsize=16, labelpad=15)
                plt.ylabel("Transit Depth (%)", fontsize=16, labelpad=15)
                axes = plt.gca()
                axes.set_xlim(wavelimlo, wavelimhi)
    
                # Overlay rebinned model spectrum and define plot filenames
                plt.plot(x1, rebinned_model*100 , 'tab:gray', lw=4)
                plt.tick_params(labelsize=14)
                plt.grid(b=True, which='major', linestyle='--')
    
                # Plot simulated transmission spectrum 
                yerr = e100
                plt.errorbar(x1, y100, yerr=yerr, fmt='o', ecolor='k',
                             mfc='k', mec='k') # showing 1-sigma errors
                plt.legend(['Model', str(Instrument)], loc=2, fontsize=16)
        
                # Save plot to .svg file
                if j == 0: 
                    plt.savefig(fname1[0], overwrite=True, bbox_inches='tight')
    
                elif j == 1:
                    plt.savefig(fname1[1], overwrite=True, bbox_inches='tight')
                    
                plt.close()    
                       
    
                ###############################################################
                # Determine the number of transits needed for strong detection
                ###############################################################    
                
                # Extract single transit simulation values with NO random noise 
                x1trans, y1trans, e1trans = single_transit['FinalSpectrum']['wave'], single_transit['FinalSpectrum']['spectrum'], single_transit['FinalSpectrum']['error_w_floor']


        
                for idx, ntr in enumerate(ints):
                    
                    # BIC loop to account for jitter in the noisy data.. 
                    for nrun in range(0, Npdxo): 
                
                        np.random.seed()
                        multi_trans_noise = e1trans/np.sqrt(ntr) 
                        
                        # add in noise floor
                        multi_trans_noise[multi_trans_noise 
                                          < nfloor/10**6] = nfloor/10**6 
                        
                        # add randomness
                        rand_noise= multi_trans_noise*(np.random.randn(len(x1trans))) 
                        
                        # this rand_noise will get smaller and smaller 
                        # as your ntr increases 
                        simulated_spectrum = y1trans + rand_noise  
                        
                        x, y, e = x1trans, simulated_spectrum, multi_trans_noise       
            
                                                   
                        x1 = np.ndarray.flatten(np.array(x))
                        y100 = np.ndarray.flatten(np.asarray(y)*100)
                        e100 = np.ndarray.flatten(np.asarray(e)*100)
    
    
                        #######################################################                                           
                        # Determine delta BIC, comparing sim data to ExoT model 
                        # and flat line
                        #######################################################
                        
                        # Model spectrum case
                        ndof = 4
                        
                        # Determine model chi-squared and reduced chi-squared
                        chi_squared = np.sum(((y100-y_rebin)/e100)**2)
                        reduced_chi_squared = chi_squared/(len(x1)-ndof)
                        
                        # Determine error factor to drive reduced 
                        # chi-squared to 1.
                        errfactor = np.sqrt(reduced_chi_squared)
                    
                        # Determine rescaled chi-squared for this model
                        sig_comb = np.sqrt((errfactor*e100)**2)
                        chi_squared = np.sum(((y100-y_rebin)/sig_comb)**2)
                        reduced_chi_squared = chi_squared/(len(x1)-ndof)
                                       
                        # Compute BIC for model spectrum case
                        BICMS = chi_squared + ndof*np.log(len(x1))
                                          
                        
                        # Flat-line spectrum case
                        ndof = 1
                        
                        # Determine chi-squared for flat-line model
                        chi_squared = np.sum(((y100
                                               -np.median(y100))/sig_comb)**2)
                                     
                        # Compute BIC for flat-line spectrum case
                        BICFL = chi_squared + ndof*np.log(len(x1))
                        
                        xnt[idx, nrun] = ntr
                        dBIC[idx, nrun] = BICFL - BICMS
                        
                    
                    # Compute mean dBIC for all transit cases    
                    dBIC_mean[idx] = np.mean(dBIC[idx])
                      
                    # Compute std dev of the dBIC for all transit cases    
                    dBIC_sdev[idx] = np.std(dBIC[idx])
                    
            
                # Define the base data for interpolation:
                for ii in range(0, np.max(idx)+1):
                    xx[ii] = xnt[ii, 0]
                yy = dBIC_mean - dBIC_sdev 
           
            
                # Force yy and dBIC_mean to increase monotonically
                for iii in range(0, np.max(idx)):
                    if yy[iii+1] < yy[iii]:
                        yy[iii+1] = yy[iii]
                        
                        
                for iii in range(0, np.max(idx)):
                    if dBIC_mean[iii+1] < dBIC_mean[iii]:
                        dBIC_mean[iii+1] = dBIC_mean[iii]  
                        
                
                # First adjust smoothing parameter for interpolation routine
                if dBIC_mean[np.max(idx)] > 175:
                    smooth = 0
                else:
                    smooth = 2
                    
                # Define interpolation function for 
                # mean dBIC (less 1 sigma) vs #transits (soft smoothing)    
                dbicvstr1 = UnivariateSpline(xx, yy, s=smooth)    


                # Check if dBIC_mean (less one sigma) for the max number 
                # of transits is < dBIC_thresh
                if (dbicvstr1(xx[np.max(idx)])) < dBIC_thresh:
                    
                    if j == 0:
                        print('               mean dBIC (less one sigma) <', 
                              dBIC_thresh, 'for', str(np.max(ints)), 
                              'transits ---> no detection')
                        print('\n') 
                        
                        # Set flags
                        nt_lo[n-Nrow_start] = 999  
                        tT_lo[n-Nrow_start] = 9990
    
    
                    elif j == 1:
                        print('               mean dBIC (less one sigma) <',
                              dBIC_thresh, 'for', str(np.max(ints)), 
                              'transits ---> no detection')
                        print('\n') 
                        
                        # Set flags
                        nt_hi[n-Nrow_start] = 999  
                        tT_hi[n-Nrow_start] = 9990   
    
    
                # Check if dBIC_mean (less one sigma) for a 
                # single transit is > dBIC_thresh
                elif (dbicvstr1(xx[0])) > dBIC_thresh:
                   
                    if j == 0:
                        nt_lo[n-Nrow_start] = 1
                        print('               min # of transits for mean'
                              ' dBIC (less one sigma) >', dBIC_thresh, ': ', 1)
                        print('\n') 
                        
                        # Determine total instrument obs cycle for detection
                        tT_lo[n-Nrow_start] = ((2 * tdur_sec/3600 
                             + 2.4)*1.16)*nt_lo[n-Nrow_start]
                        
                    
                    elif j == 1:
                        nt_hi[n-Nrow_start] = 1
                        print('               min # of transits for mean'
                              ' dBIC (less one sigma) >', dBIC_thresh, ': ', 1)
                        print('\n') 
                        
                        # Determine total instrument obs cycle for detection
                        tT_hi[n-Nrow_start] = ((2 * tdur_sec/3600 
                             + 2.4)*1.16)*nt_hi[n-Nrow_start]
                                   
                
                
                ###############################################################
                # Interpolate dBIC vs #transits 
                ###############################################################
                else: 
                                            
                    # Set up plot axes, etc.
                    fignum = 2
                    plt.figure(fignum, figsize = (17, 9))
                    axes = plt.gca()
                    plt.xlabel('# transits', fontsize=16, labelpad=15)
                    plt.ylabel("dBIC", fontsize=16, labelpad=15)
                    
                    # Plot dBIC threshold
                    ythresh = np.linspace(dBIC_thresh, dBIC_thresh, len(ints)) 
                    plt.plot(xx, ythresh, 'k--', lw=3, color='tab:gray')
                    
                    # Plot dBIC_mean vs #transits, grid pts.
                    plt.plot(xx, dBIC_mean, 'k.')
                    
                    # Define interpolation function for 
                    # mean dBIC vs #transits (soft smoothing)
                    dbicvstr = UnivariateSpline(xx, dBIC_mean, s=smooth)
                    
                    # Define finer grid for plot interpolation
                    xxfine = np.linspace(1, np.max(ints), 5*np.max(ints))

                    # Plot dBIC_mean vs #transits interp.
                    plt.plot(xxfine, dbicvstr(xxfine), 'k:', lw=1)
                                                                                  
                    # Plot mean dBIC (less one sigma) vs #transits interp.
                    plt.plot(xxfine, dbicvstr1(xxfine), 'k-', lw=2)
                    
                    
                    plt.tick_params(labelsize=14)
                    plt.grid(b=True, which='major', linestyle='-.')
                    plt.legend(['very strong detection threshold', 
                                r'mean dBIC grid pts.', 'mean dBIC',
                                r'mean dBIC (less 1$\sigma$)'], 
                                 loc=4, fontsize=16)
    
                    fname2 = [os.getcwd() + '/dBIC/' + 'dBIC_'
                              + str(n) + '_lo_metal.svg',
                              os.getcwd() + '/dBIC/' + 'dBIC_'
                              + str(n) + '_hi_metal.svg']
                    
                    # Save plot to .svg file
                    if j == 0:    
                        plt.savefig(fname2[0], overwrite=True, 
                                    bbox_inches='tight')
                        
                    elif j == 1:
                        plt.savefig(fname2[1], overwrite=True, 
                                    bbox_inches='tight')
        
                    plt.close()
                              
                        
                    ###########################################################
                    # Using zero-finder, find min # transits 
                    # for dBIC > dBIC_thresh
                    ###########################################################  
                                                  
                    def froot(x):
                        return dbicvstr1(x) - dBIC_thresh
                    
                    
                    root = optimize.brentq(froot, a=1.0, b=np.max(ints), 
                                           maxiter=5000)
                    
                    if j == 0:
                        nt_lo[n-Nrow_start] = np.ceil(root)
                        print('                 min # of transits for mean'
                              ' dBIC (less one sigma) >', dBIC_thresh,
                              ': ', nt_lo[n-Nrow_start])
                        print('\n') 
                        
                        # Determine total instrument obs cycle for detection
                        tT_lo[n-Nrow_start] = ((2 * tdur_sec/3600 
                             + 2.4)*1.16)*nt_lo[n-Nrow_start]
                        
                        # Check if number of transits are observable in 10 yrs
                        if pdxinput[ntransfer, 9] < nt_lo[n-Nrow_start]:
                            print('               # of transits needed for'
                                  ' detection exceed those available'
                                  ' in 10yr mission')
                            print('\n')
                            tT_lo[n-Nrow_start] = 9970
                                                       
                    elif j == 1:
                        nt_hi[n-Nrow_start] = np.ceil(root)
                        print('               min # of transits for mean'
                              ' dBIC (less one sigma) >', dBIC_thresh,
                              ': ', nt_hi[n-Nrow_start])
                        print('\n') 
                                                
                        # Determine total instrument obs cycle for detection
                        tT_hi[n-Nrow_start] = ((2 * tdur_sec/3600 + 
                             2.4)*1.16)*nt_hi[n-Nrow_start]
                        
                        # Check if number of transits are observable in 10 yrs
                        if pdxinput[ntransfer, 9] < nt_hi[n-Nrow_start]:
                            print('               # of transits needed for'
                                  ' detection exceed those observable'
                                  ' in a 10yr mission')
                            print('\n')
                            tT_hi[n-Nrow_start] = 9970
                                           
                        
        # Print warning note if Jmag is too low     
        else:
            print('           Warning: Jmag too low for this target,'
                  ' detector saturation \n')
            nt_lo[n-Nrow_start] = 998     # if Jmag < Jmagsat for NIRSpec G235M
            tT_lo[n-Nrow_start] = 9980
            nt_hi[n-Nrow_start] = 998     # if Jmag < Jmagsat for NIRSpec G235M
            tT_hi[n-Nrow_start] = 9980
            
            
    # Print warning note if Rp is too high     
    else:
        print('           Warning: This target has Rp > 10;'
              ' unreliable mass est. \n')
        nt_lo[n-Nrow_start] = 996         # if Rp > 10 Rearth
        tT_lo[n-Nrow_start] = 9960
        nt_hi[n-Nrow_start] = 996         # if Rp > 10 Rearth
        tT_hi[n-Nrow_start] = 9960        
            
            
           
###############################################################################
# Consolidate summary data
###############################################################################

# Initialize summary_data array 
Nrows = Nrow_end - Nrow_start + 1                 
summary_data = np.zeros((Nrows, 26))                               

for i in range(0, Nrow_end-Nrow_start+1):                                                    
    
    n = i + Nrow_start                                                      
    ntransfer = n - int(pdxinput[0, 0])                                

    summary_data[i, 0]  = pdxinput[ntransfer, 10]                # Cat        
    summary_data[i, 1]  = 1                                      # Rank      
    summary_data[i, 2]  = n                                      # Target     
    summary_data[i, 3]  = survey_data[i+Nrow_start-1, 0]         # RAdeg     
    summary_data[i, 4]  = survey_data[i+Nrow_start-1, 1]         # DEdeg      
    summary_data[i, 5]  = survey_data[i+Nrow_start-1, 2]         # Rp
    summary_data[i, 6]  = survey_data[i+Nrow_start-1, 3]         # Per
    summary_data[i, 7]  = survey_data[i+Nrow_start-1, 4]         # S
    summary_data[i, 8]  = survey_data[i+Nrow_start-1, 5]         # K
    summary_data[i, 9]  = survey_data[i+Nrow_start-1, 6]         # R*
    summary_data[i, 10] = survey_data[i+Nrow_start-1, 7]         # Teff
    summary_data[i, 11] = survey_data[i+Nrow_start-1, 10]        # Jmag
    summary_data[i, 12] = pdxinput[ntransfer, 1]                 # a          
    summary_data[i, 13] = pdxinput[ntransfer, 2]                 # Teq
    summary_data[i, 14] = pdxinput[ntransfer, 6]                 # gs
    summary_data[i, 15] = pdxinput[ntransfer, 3]                 # Mmed
    summary_data[i, 18] = pdxinput[ntransfer, 7]                 # Rhomed
    summary_data[i, 19] = pdxinput[ntransfer, 8]                 # tdur
    summary_data[i, 20] = pdxinput[ntransfer, 9]                 # nt10yr
    summary_data[i, 21] = nt_lo[n-Nrow_start]                    # nt_lo      
    summary_data[i, 22] = nt_hi[n-Nrow_start]                    # nt_hi
    summary_data[i, 23] = tT_lo[n-Nrow_start]                    # tT_lo      
    summary_data[i, 24] = tT_hi[n-Nrow_start]                    # tT_hi      
    summary_data[i, 25] = (tT_lo[n-Nrow_start] + 
                           tT_hi[n-Nrow_start])/2                # tT_avg      



###############################################################################
# Write transfer_data to Pdxo output file
###############################################################################
summary_data_old = np.loadtxt('Pdxo_Output.txt', ndmin=2, skiprows=0)
targets_old = np.max(summary_data_old, axis=0)
last_target_old = targets_old[2]

targets_new = np.min(summary_data, axis=0)
first_target_new = targets_new[2]

if first_target_new == last_target_old + 1:
    with open('Pdxo_Output.txt', 'ba') as f: 
        np.savetxt(f, summary_data)
else:
    np.savetxt('Pdxo_Output.txt', summary_data)


raise SystemExit