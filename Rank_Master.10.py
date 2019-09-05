

"""
@author: Charles Fortenbach

Code associated with manuscript submitted to PASP (TBD)

Title: A FRAMEWORK FOR OPTIMIZING EXOPLANET TARGET SELECTION FOR
THE JAMES WEBB SPACE TELESCOPE

Subprogram:  Rank_Master.py

Version:  1.1

This project is licensed under the GNU GPLv3 License.

"""

import numpy as np
import textwrap
import datetime
import time
import os



###############################################################################
# Read in JET input file
###############################################################################
JET_params = np.genfromtxt('JET_Input.txt', dtype=str, delimiter=', ', 
             skip_header=6)

# Target catalog starting and ending rows:
Nrow_start = int(JET_params[0, 1])
Nrow_end = int(JET_params[1, 1])


# JWST instrument/mode
Instrument = JET_params[2, 1]

# Wevelength range for this instrument (microns):
wavelimlo = float(JET_params[3, 1])   
wavelimhi = float(JET_params[4, 1])

# Brightness limit (Jmag) at 100% FW as f(Teff):   
Jmag_10000K_100FW = float(JET_params[5, 1])
Jmag_5000K_100FW = float(JET_params[6, 1])
Jmag_2500K_100FW = float(JET_params[7, 1])

# Detector linearity limit as fraction of FW (%):
percentFW = float(JET_params[8, 1])

# Instrument multi-transit residual noise floor (ppm):
nfloor = float(JET_params[9, 1])

# R value of simulation:
Res = int(JET_params[10, 1])

# Number of samples for each dBIC vs ntr g pt. (Npdxo):
Npdxo = int(JET_params[11, 1])

# Detection threshold (dBIC):
dBIC_thresh = int(JET_params[12, 1])


# Eq. of State (lo metal atm):
EOS_lo = JET_params[13, 1]

# Cloud_lo (Pa):
Cloud_lo = float(JET_params[14, 1])

# Eq. of State  (hi metal atm):
EOS_hi = JET_params[15, 1]

# Cloud_hi (Pa):
Cloud_hi = float(JET_params[16, 1])


# Observing time elements incl. overheads:
# Out of transit factor (% tdur):      
Out_trans_factor = float(JET_params[17, 1])

# + Out of transit "timing tax"(sec): 
Timing_tax = float(JET_params[18, 1])

# Slew duration avg. (sec):
Slew_time = float(JET_params[19, 1])

# SAMs: small angle maneuvers (sec):
SAMs = float(JET_params[20, 1])

# GS Acq: guide star acquisition(s)(sec):
GS_acq = float(JET_params[21, 1])

# Targ Acq: target acquisition if any (sec):
Targ_acq = float(JET_params[22, 1])

# Exposure Ovhd: factor 1: 
Exp_Ovhd_factor1 = float(JET_params[23, 1])

# Exposure Ovhd: factor 2 (sec): 
Exp_Ovhd_factor2 = float(JET_params[24, 1])

# Mech: mechanism movements (sec):
Mech = float(JET_params[25, 1])

# OSS: Onboard Script System compilation (sec):
OSS = float(JET_params[26, 1])

# MSA: NIRSpec MSA configuration (sec):  
MSA = float(JET_params[27, 1])

# IRS2: NIRSpec IRS2 Detector Mode setup (sec):
IRS2 = float(JET_params[28, 1])

# Visit Ovhd: visit cleanup activities (sec): 
Visit_Ovhd = float(JET_params[29, 1])

# Obs Ovhd factor (%):
Obs_Ovhd_factor = float(JET_params[30, 1])

# DS Ovhd (sec):
DS_Ovhd = float(JET_params[31, 1])


# Set Run flags
RunExoT = JET_params[32, 1]
RunPdxo = JET_params[33, 1]
RunRank = JET_params[34, 1]

   
        
###############################################################################
# Check Run flag
###############################################################################
if RunRank == 'N':
    raise SystemExit
  
    
    
###############################################################################
# Print notification of running Pdxo_Master
###############################################################################
print('\n')                      
print(" Running Rank_Master:")



###############################################################################
# Read in survey data file 
###############################################################################
summary_data = np.loadtxt('Pdxo_Output.txt', ndmin=2, skiprows=0)


# Reset target indices
Nrow_start = int(summary_data[0, 2])
targets = np.max(summary_data, axis=0)
Nrow_end = int(targets[2])
    


###############################################################################
# Initialize arrays
###############################################################################
start = np.zeros(7, dtype=int)
end = np.zeros(7, dtype=int)



###############################################################################
# Write summary table for un-ranked targets to txt file in working directory
###############################################################################
now = datetime.datetime.now()
timestr = time.strftime("%Y%m%d-%H%M%S")

fout = open(os.getcwd() + '/JET_Smry_tables/' +'JET_Smry_'
            + timestr + '_un-ranked.txt', "w", newline='')

fout.write('JET run ID:  ' + str(now) + '\n' + '\n')

fout.write('Input Summary: ' + '\n' + '\n')
fout.write(' Parameter                                          ' 
           + 'Input Value' + '\n')
fout.write(' ---------                                          ' 
           + '-----------' + '\n')
fout.write(' Starting target for this sequence of runs:         ' 
           + str(Nrow_start) + '\n')
fout.write(' Ending target from catalog:                        ' 
           + str(Nrow_end) + '\n' + '\n')

fout.write(' JWST instrument:                                   ' 
           + Instrument + '\n')
fout.write(' Wavelength short limit (microns):                  ' 
           + str(wavelimlo) + '\n')
fout.write(' Wavelength long limit (microns):                   ' 
           + str(wavelimhi) + '\n')
fout.write(' Jmag limit (Teff = 10000K):                        ' 
           + str(Jmag_10000K_100FW) + '\n')
fout.write(' Jmag limit (Teff = 5000K):                         ' 
           + str(Jmag_5000K_100FW) + '\n')
fout.write(' Jmag limit (Teff = 2500K):                         ' 
           + str(Jmag_2500K_100FW) + '\n')
fout.write(' Detector linear response limit (% FW):             ' 
           + str(percentFW) + '\n')
fout.write(' Noise floor (nfloor) (ppm):                        ' 
           + str(nfloor) + '\n')
fout.write(' R value of sim (Res):                              ' 
           + str(Res) + '\n')
fout.write(' Number of samples for dBIC vs ntr grid pts.:       ' 
           + str(Npdxo) + '\n')
fout.write(' Detection threshold (dBIC):                        ' 
           + str(dBIC_thresh) + '\n' + '\n')

fout.write(' Eq. of State (lo metal atm):                       ' 
           + EOS_lo + '\n')
fout.write(' Cloud_lo (Pa):                                     ' 
           + str(Cloud_lo) + '\n')
fout.write(' Eq. of State (hi metal atm):                       ' 
           + EOS_hi + '\n')
fout.write(' Cloud_hi (Pa):                                     ' 
           + str(Cloud_hi) + '\n' + '\n')

fout.write(' Out of transit factor (% tdur):                    ' 
           + str(Out_trans_factor) + '\n')
fout.write(' + Out of transit "timing tax" (sec):               ' 
           + str(Timing_tax) + '\n')
fout.write(' Slew duration avg. (sec):                          ' 
           + str(Slew_time) + '\n')
fout.write(' SAMs: small angle maneuvers (sec):                 ' 
           + str(SAMs) + '\n')
fout.write(' GS Acq: guide star acquisition(s)(sec):            ' 
           + str(GS_acq) + '\n')
fout.write(' Targ Acq: target acquisition if any (sec):         ' 
           + str(Targ_acq) + '\n')
fout.write(' Exposure Ovhd: factor 1:                           ' 
           + str(Exp_Ovhd_factor1) + '\n')
fout.write(' Exposure Ovhd: factor 2 (sec):                     ' 
           + str(Exp_Ovhd_factor2) + '\n')
fout.write(' Mech: mechanism movements (sec):                   ' 
           + str(Mech) + '\n')
fout.write(' OSS: Onboard Script System compilation (sec):      ' 
           + str(OSS) + '\n')
fout.write(' MSA: NIRSpec MSA configuration (sec):              ' 
           + str(MSA) + '\n')
fout.write(' IRS2: NIRSpec IRS2 Detector Mode setup (sec):      ' 
           + str(IRS2) + '\n')
fout.write(' Visit Ovhd: visit cleanup activities (sec):        ' 
           + str(Visit_Ovhd) + '\n')
fout.write(' Obs Ovhd factor (%):                               ' 
           + str(Obs_Ovhd_factor) + '\n')
fout.write(' DS_Ovhd (sec):                                     ' 
           + str(DS_Ovhd) + '\n' + '\n')

fout.write(textwrap.dedent("""\

 Output Summary (un-ranked):
                      
 Col  Label     Explanations
 --------------------------------------------------------------------------
 ---  ranking:
   0  Cat       Planet demographic category (1 - 7)
  
 ---  survey data:
   1  Targ      Survey table index (starting with 1)
   2  RAdeg     Right Ascension in decimal degrees (J2000)
   3  DEdeg     Declination in decimal degrees (J2000)
   4  Rp        Planetary radius (mean) in Earth units (R/Rearth)
   5  Per       Period (days)
   6  S         Planetary insolation in Earth units (S/Searth)
   7  R*        Stellar radius (Solar radii)
   8  Teff      Stellar effective temperature (K)
   9  Jmag      Apparent J band magnitude (mag)
  
 ---  computed values:
  10  a         Semi-major axis of planet orbit (Solar radii)
  11  Teq       Planet equilibrium temperature (K)
  12  gs        Planet surface gravity (m/s^2)
  13  Mmed      Planet mass (median) in Earth units (M/Mearth)
  14  tdur      Transit duration (hrs)
  15  nt10yr    Number of transits observable in 10yr mission
  16  nt_lo     Number of transits needed for detection; lo-metal atm
  17  nt_hi     Number of transits needed for detection; hi-metal atm
  18  tT_lo     Total observ. cycle (in/out of transit); lo-metal atm (hrs)
  19  tT_hi     Total observ. cycle (in/out of transit); hi-metal atm (hrs)
  21  tT_avg    Total observ. cycle for 'average' atm detection (hrs)
 """))
    
fout.write('\n')
fout.write('\n')


# Define headers and formats for output
head1 = ("               |----------------------------------survey data-------"
         "-------------------------|"
         "---------------------------------------------computed values--------"
         "-----------------------------|\n")
head2 = ("  Cat            Targ    RAdeg     DEdeg      Rp      Per          "
         "S      R*   Teff     Jmag          a    Teq        gs      "
         "Mmed     tdur  nt10yr  nt_lo  nt_hi     "
         "tT_lo     tT_hi    tT_avg\n")
fmt1 = (" {Cat:4.0f} {n:15.0f} {RAdeg:8.4f} {DEdeg:9.4f} {Rp:7.4f} {Per:8.4f}"
        " {S:10.3f} {Rstar:7.4f} {Teff:6.0f} {Jmag:8.3f} {a:10.4f} {Teq:6.0f}"
        " {gs:9.4f} {Mmed:9.4f} {tdur:8.4f} {nt10yr:7.0f} {nt_lo:6.0f}"
        " {nt_hi:6.0f} {tT_lo:9.1f} {tT_hi:9.1f} {tT_avg:9.1f}\n")

fout.write(head1)
fout.write(head2)

for i in range(0, Nrow_end - Nrow_start + 1):           
               
    fout.write(fmt1.format(
         Cat      = summary_data[i, 0],
         Rank     = summary_data[i, 1], 
         n        = summary_data[i, 2], 
         RAdeg    = summary_data[i, 3], 
         DEdeg    = summary_data[i, 4], 
         Rp       = summary_data[i, 5], 
         Per      = summary_data[i, 6], 
         S        = summary_data[i, 7], 
         K        = summary_data[i, 8],
         Rstar    = summary_data[i, 9], 
         Teff     = summary_data[i, 10], 
         Jmag     = summary_data[i, 11], 
         a        = summary_data[i, 12],
         Teq      = summary_data[i, 13],                   
         gs       = summary_data[i, 14], 
         Mmed     = summary_data[i, 15], 
         Rhomed   = summary_data[i, 18], 
         tdur     = summary_data[i, 19], 
         nt10yr   = summary_data[i, 20],
         nt_lo    = summary_data[i, 21],
         nt_hi    = summary_data[i, 22],
         tT_lo    = summary_data[i, 23],
         tT_hi    = summary_data[i, 24],
         tT_avg   = summary_data[i, 25]))

    
fout.write(textwrap.dedent("""\


                           
Notes: 
    "995 and 9950.0"   tT_avg undetermined since hi_metal atm not detected 
    "996 and 9960.0"   Rp > 10:  Mmed estimates are unreliable
    "997 and 9970.0"   number of transits needed for detection exceeds"""
    + """ those observable in 10 yr mission
    "998 and 9980.0"   target Jmag is below (brighter than) the saturation"""
    + """ threshold of this instrument/mode
    "999 and 9990.0"   dBIC is < dBIC threshold for this target: """  
    + """no detection for this atm assumption    
    
 """))    

fout.close()


###############################################################################
# Sort summary data by category and then total observation time
###############################################################################
print('\n') 
print(" Sorting and ranking targets . . .\n\n" ) 

summary_data = sorted(summary_data, key=lambda x: (x[0], x[25]))
summary_data = np.asarray(summary_data)

for i in range(0, Nrow_end - Nrow_start):
    if summary_data[i+1, 0] == summary_data[i, 0]:     # if Cat is the same     
        summary_data[i+1, 1] = summary_data[i, 1] + 1  # then increment rank
    else:
        summary_data[i+1, 1] = 1   # otherwise reset the rank to one, and loop
       
        

###############################################################################
# Write summary table to txt file in working directory
###############################################################################
fout = open(os.getcwd() + '/JET_Smry_tables/' +'JET_Smry_' 
            + timestr + '.txt', "w", newline='')

fout.write('JET run ID:  ' + str(now) + '\n' + '\n')

fout.write('Input Summary: ' + '\n' + '\n')
fout.write(' Parameter                                          ' 
           + 'Input Value' + '\n')
fout.write(' ---------                                          ' 
           + '-----------' + '\n')
fout.write(' Starting target for this sequence of runs:         ' 
           + str(Nrow_start) + '\n')
fout.write(' Ending target from catalog:                        ' 
           + str(Nrow_end) + '\n' + '\n')

fout.write(' JWST instrument:                                   ' 
           + Instrument + '\n')
fout.write(' Wavelength short limit (microns):                  ' 
           + str(wavelimlo) + '\n')
fout.write(' Wavelength long limit (microns):                   ' 
           + str(wavelimhi) + '\n')
fout.write(' Jmag limit (Teff = 10000K):                        ' 
           + str(Jmag_10000K_100FW) + '\n')
fout.write(' Jmag limit (Teff = 5000K):                         ' 
           + str(Jmag_5000K_100FW) + '\n')
fout.write(' Jmag limit (Teff = 2500K):                         ' 
           + str(Jmag_2500K_100FW) + '\n')
fout.write(' Detector linear response limit (% FW):             ' 
           + str(percentFW) + '\n')
fout.write(' Noise floor (nfloor) (ppm):                        ' 
           + str(nfloor) + '\n')
fout.write(' R value of sim (Res):                              ' 
           + str(Res) + '\n')
fout.write(' Number of samples for dBIC vs ntr grid pts.:       ' 
           + str(Npdxo) + '\n')
fout.write(' Detection threshold (dBIC):                        ' 
           + str(dBIC_thresh) + '\n' + '\n')

fout.write(' Eq. of State (lo metal atm):                       ' 
           + EOS_lo + '\n')
fout.write(' Cloud_lo (Pa):                                     ' 
           + str(Cloud_lo) + '\n')
fout.write(' Eq. of State (hi metal atm):                       ' 
           + EOS_hi + '\n')
fout.write(' Cloud_hi (Pa):                                     ' 
           + str(Cloud_hi) + '\n' + '\n')

fout.write(' Out of transit factor (% tdur):                    ' 
           + str(Out_trans_factor) + '\n')
fout.write(' + Out of transit "timing tax" (sec):               ' 
           + str(Timing_tax) + '\n')
fout.write(' Slew duration avg. (sec):                          ' 
           + str(Slew_time) + '\n')
fout.write(' SAMs: small angle maneuvers (sec):                 ' 
           + str(SAMs) + '\n')
fout.write(' GS Acq: guide star acquisition(s)(sec):            ' 
           + str(GS_acq) + '\n')
fout.write(' Targ Acq: target acquisition if any (sec):         ' 
           + str(Targ_acq) + '\n')
fout.write(' Exposure Ovhd: factor 1:                           ' 
           + str(Exp_Ovhd_factor1) + '\n')
fout.write(' Exposure Ovhd: factor 2 (sec):                     ' 
           + str(Exp_Ovhd_factor2) + '\n')
fout.write(' Mech: mechanism movements (sec):                   ' 
           + str(Mech) + '\n')
fout.write(' OSS: Onboard Script System compilation (sec):      ' 
           + str(OSS) + '\n')
fout.write(' MSA: NIRSpec MSA configuration (sec):              ' 
           + str(MSA) + '\n')
fout.write(' IRS2: NIRSpec IRS2 Detector Mode setup (sec):      ' 
           + str(IRS2) + '\n')
fout.write(' Visit Ovhd: visit cleanup activities (sec):        ' 
           + str(Visit_Ovhd) + '\n')
fout.write(' Obs Ovhd factor (%):                               ' 
           + str(Obs_Ovhd_factor) + '\n')
fout.write(' DS_Ovhd (sec):                                     ' 
           + str(DS_Ovhd) + '\n' + '\n')

fout.write(textwrap.dedent("""\

 Output Summary (fully sorted and ranked):
     
 Col  Label     Explanations
 --------------------------------------------------------------------------
      Category  Planet demographic category (1 - 7)
      
 ---  ranking:      
   0  Rank      Rank within Cat based on total observ. cycle for detection
  
 ---  survey data:
   1  Targ      Survey table index (starting with 1)
   2  RAdeg     Right Ascension in decimal degrees (J2000)
   3  DEdeg     Declination in decimal degrees (J2000)
   4  Rp        Planetary radius (mean) in Earth units (R/Rearth)
   5  Per       Period (days)
   6  S         Planetary insolation in Earth units (S/Searth)
   7  R*        Stellar radius (Solar radii)
   8  Teff      Stellar effective temperature (K)
   9  Jmag      Apparent J band magnitude (mag)
  
 ---  computed values:
  10  a         Semi-major axis of planet orbit (Solar radii)
  11  Teq       Planet equilibrium temperature (K)
  12  gs        Planet surface gravity (m/s^2)
  13  Mmed      Planet mass (median) in Earth units (M/Mearth)
  14  tdur      Transit duration (hrs)
  15  nt10yr    Number of transits observable in 10yr mission
  16  nt_lo     Number of transits needed for detection; lo-metal atm
  17  nt_hi     Number of transits needed for detection; hi-metal atm
  18  tT_lo     Total observ. cycle (in/out of transit); lo-metal atm (hrs)
  19  tT_hi     Total observ. cycle (in/out of transit); hi-metal atm (hrs)
  20  tT_avg    Total observ. cycle for 'average' atm detection (hrs)
 """))

    
# Define headers and formats for output
head1 = ("      |----------------------------------survey data----------------"
         "----------------|"
         "--------------------------------------------computed values---------"
         "-----------------------------|\n")
head2 = (" Rank   Targ    RAdeg     DEdeg      Rp      Per          S      R* "
         "  Teff     Jmag          a    Teq        gs      "
         "Mmed     tdur  nt10yr  nt_lo  nt_hi     "
         "tT_lo     tT_hi    tT_avg\n")
fmt1  = (" {Rank:4.0f} {n:6.0f} {RAdeg:8.4f} {DEdeg:9.4f} {Rp:7.4f} {Per:8.4f}"
        " {S:10.3f} {Rstar:7.4f} {Teff:6.0f} {Jmag:8.3f} {a:10.4f} {Teq:6.0f}"
        " {gs:9.4f} {Mmed:9.4f} {tdur:8.4f} {nt10yr:7.0f} {nt_lo:6.0f}"
        " {nt_hi:6.0f} {tT_lo:9.1f} {tT_hi:9.1f} {tT_avg:9.1f}\n")
head3 = (" Rank   Targ    RAdeg     DEdeg      Rp      Per          S      R* "
         "  Teff     Jmag          a    Teq        gs      "
         "Mmed     tdur  nt10yr\n")
fmt2  = (" Min:                           {Rp:7.4f} {Per:8.4f}"
        " {S:10.3f} {Rstar:7.4f} {Teff:6.0f} {Jmag:8.3f} {a:10.4f} {Teq:6.0f}"
        " {gs:9.4f} {Mmed:9.4f} {tdur:8.4f} {nt10yr:7.0f}\n")
fmt3  = (" Med:                           {Rp:7.4f} {Per:8.4f}"
        " {S:10.3f} {Rstar:7.4f} {Teff:6.0f} {Jmag:8.3f} {a:10.4f} {Teq:6.0f}"
        " {gs:9.4f} {Mmed:9.4f} {tdur:8.4f} {nt10yr:7.0f}\n")
fmt4  = ("Mean:                           {Rp:7.4f} {Per:8.4f}"
        " {S:10.3f} {Rstar:7.4f} {Teff:6.0f} {Jmag:8.3f} {a:10.4f} {Teq:6.0f}"
        " {gs:9.4f} {Mmed:9.4f} {tdur:8.4f} {nt10yr:7.0f}\n")
fmt5  = (" Max:                           {Rp:7.4f} {Per:8.4f}"
        " {S:10.3f} {Rstar:7.4f} {Teff:6.0f} {Jmag:8.3f} {a:10.4f} {Teq:6.0f}"
        " {gs:9.4f} {Mmed:9.4f} {tdur:8.4f} {nt10yr:7.0f}\n")

head5 = "Full Target Set:\n" 
head4 = ("        Targ    RAdeg     DEdeg      Rp      Per          S      R* "
         "  Teff     Jmag          a    Teq        gs      "
         "Mmed     tdur  nt10yr\n")
fmt6  = (" Min: {n:6.0f} {RAdeg:8.4f} {DEdeg:9.4f} {Rp:7.4f} {Per:8.4f}"
        " {S:10.3f} {Rstar:7.4f} {Teff:6.0f} {Jmag:8.3f} {a:10.4f} {Teq:6.0f}"
        " {gs:9.4f} {Mmed:9.4f} {tdur:8.4f} {nt10yr:7.0f}\n")
fmt7  = (" Med:                           {Rp:7.4f} {Per:8.4f}"
        " {S:10.3f} {Rstar:7.4f} {Teff:6.0f} {Jmag:8.3f} {a:10.4f} {Teq:6.0f}"
        " {gs:9.4f} {Mmed:9.4f} {tdur:8.4f} {nt10yr:7.0f}\n")
fmt8  = ("Mean:                           {Rp:7.4f} {Per:8.4f}"
        " {S:10.3f} {Rstar:7.4f} {Teff:6.0f} {Jmag:8.3f} {a:10.4f} {Teq:6.0f}"
        " {gs:9.4f} {Mmed:9.4f} {tdur:8.4f} {nt10yr:7.0f}\n")
fmt9  = (" Max: {n:6.0f} {RAdeg:8.4f} {DEdeg:9.4f} {Rp:7.4f} {Per:8.4f}"
        " {S:10.3f} {Rstar:7.4f} {Teff:6.0f} {Jmag:8.3f} {a:10.4f} {Teq:6.0f}"
        " {gs:9.4f} {Mmed:9.4f} {tdur:8.4f} {nt10yr:7.0f}\n")


# Set start and end rows for each category
for j in range(0, 7): 
    for i in range(0, Nrow_end - Nrow_start + 1):        
        if np.int(summary_data[i, 0]) == j + 1:            
            if np.int(summary_data[i, 1]) == 1:               
                start[j] = i
                    
for j in range(0, 7): 
    for i in range(0, Nrow_end - Nrow_start):        
        if np.int(summary_data[i, 0]) == j + 1:            
            if np.int(summary_data[i+1, 1]) <= np.int(summary_data[i, 1]):               
                end[j] = i

for j in range(0, 7): 
    for i in range(0, Nrow_end - Nrow_start + 1): 
        if j == np.max(summary_data[:, 0])-1 and i == (Nrow_end - Nrow_start):
            end[j] = i

# Set flag for targets with tT_avg undetermined since hi_metal atm not detected
for i in range(0, Nrow_end - Nrow_start + 1):
    if summary_data[i, 25] > 4000 and summary_data[i, 25] < 6000:
        summary_data[i, 25] = 9950    


# Write summary table
for j in range(0, 7):
    fout.write('\n')
    fout.write('\n')
    fout.write('Category: ' + str(j+1)) 
    fout.write('\n')
    fout.write(head1)
    fout.write(head2)
                
    for i in range(0, Nrow_end - Nrow_start + 1):           
        
        if np.int(summary_data[i, 0]) == j + 1:
  
            fout.write(fmt1.format(
                 Cat      = summary_data[i, 0],
                 Rank     = summary_data[i, 1], 
                 n        = summary_data[i, 2], 
                 RAdeg    = summary_data[i, 3], 
                 DEdeg    = summary_data[i, 4], 
                 Rp       = summary_data[i, 5], 
                 Per      = summary_data[i, 6], 
                 S        = summary_data[i, 7], 
                 K        = summary_data[i, 8],
                 Rstar    = summary_data[i, 9], 
                 Teff     = summary_data[i, 10], 
                 Jmag     = summary_data[i, 11], 
                 a        = summary_data[i, 12],
                 Teq      = summary_data[i, 13],                   
                 gs       = summary_data[i, 14], 
                 Mmed     = summary_data[i, 15], 
                 Rhomed   = summary_data[i, 18], 
                 tdur     = summary_data[i, 19], 
                 nt10yr   = summary_data[i, 20],
                 nt_lo    = summary_data[i, 21],
                 nt_hi    = summary_data[i, 22],
                 tT_lo    = summary_data[i, 23],
                 tT_hi    = summary_data[i, 24],
                 tT_avg   = summary_data[i, 25]))
                             
   
    for i in range(0, Nrow_end - Nrow_start + 1):
    
        if np.int(summary_data[i, 0]) == j + 1:

            # Compute stats for each Category
            fout.write('\n')           
            fout.write(fmt2.format(
                 Cat      = np.min(summary_data[start[j]:end[j]+1, 0]),
                 Rank     = np.min(summary_data[start[j]:end[j]+1, 1]), 
                 n        = np.min(summary_data[start[j]:end[j]+1, 2]), 
                 RAdeg    = np.min(summary_data[start[j]:end[j]+1, 3]), 
                 DEdeg    = np.min(summary_data[start[j]:end[j]+1, 4]), 
                 Rp       = np.min(summary_data[start[j]:end[j]+1, 5]), 
                 Per      = np.min(summary_data[start[j]:end[j]+1, 6]), 
                 S        = np.min(summary_data[start[j]:end[j]+1, 7]), 
                 K        = np.min(summary_data[start[j]:end[j]+1, 8]),
                 Rstar    = np.min(summary_data[start[j]:end[j]+1, 9]), 
                 Teff     = np.min(summary_data[start[j]:end[j]+1, 10]), 
                 Jmag     = np.min(summary_data[start[j]:end[j]+1, 11]), 
                 a        = np.min(summary_data[start[j]:end[j]+1, 12]),
                 Teq      = np.min(summary_data[start[j]:end[j]+1, 13]),                   
                 gs       = np.min(summary_data[start[j]:end[j]+1, 14]), 
                 Mmed     = np.min(summary_data[start[j]:end[j]+1, 15]), 
                 Rhomed   = np.min(summary_data[start[j]:end[j]+1, 18]), 
                 tdur     = np.min(summary_data[start[j]:end[j]+1, 19]), 
                 nt10yr   = np.min(summary_data[start[j]:end[j]+1, 20]),
                 nt_lo    = np.min(summary_data[start[j]:end[j]+1, 21]),
                 nt_hi    = np.min(summary_data[start[j]:end[j]+1, 22]),
                 tT_lo    = np.min(summary_data[start[j]:end[j]+1, 23]),
                 tT_hi    = np.min(summary_data[start[j]:end[j]+1, 24]),
                 tT_avg   = np.min(summary_data[start[j]:end[j]+1, 25])))
        
            fout.write(fmt3.format(
                 Cat      = np.median(summary_data[start[j]:end[j]+1, 0]),
                 Rank     = np.median(summary_data[start[j]:end[j]+1, 1]), 
                 n        = np.median(summary_data[start[j]:end[j]+1, 2]), 
                 RAdeg    = np.median(summary_data[start[j]:end[j]+1, 3]), 
                 DEdeg    = np.median(summary_data[start[j]:end[j]+1, 4]), 
                 Rp       = np.median(summary_data[start[j]:end[j]+1, 5]), 
                 Per      = np.median(summary_data[start[j]:end[j]+1, 6]), 
                 S        = np.median(summary_data[start[j]:end[j]+1, 7]), 
                 K        = np.median(summary_data[start[j]:end[j]+1, 8]),
                 Rstar    = np.median(summary_data[start[j]:end[j]+1, 9]), 
                 Teff     = np.median(summary_data[start[j]:end[j]+1, 10]), 
                 Jmag     = np.median(summary_data[start[j]:end[j]+1, 11]), 
                 a        = np.median(summary_data[start[j]:end[j]+1, 12]),
                 Teq      = np.median(summary_data[start[j]:end[j]+1, 13]),                   
                 gs       = np.median(summary_data[start[j]:end[j]+1, 14]), 
                 Mmed     = np.median(summary_data[start[j]:end[j]+1, 15]), 
                 Rhomed   = np.median(summary_data[start[j]:end[j]+1, 18]), 
                 tdur     = np.median(summary_data[start[j]:end[j]+1, 19]), 
                 nt10yr   = np.median(summary_data[start[j]:end[j]+1, 20]),
                 nt_lo    = np.median(summary_data[start[j]:end[j]+1, 21]),
                 nt_hi    = np.median(summary_data[start[j]:end[j]+1, 22]),
                 tT_lo    = np.median(summary_data[start[j]:end[j]+1, 23]),
                 tT_hi    = np.median(summary_data[start[j]:end[j]+1, 24]),
                 tT_avg   = np.median(summary_data[start[j]:end[j]+1, 25])))
    
            fout.write(fmt4.format(
                 Cat      = np.mean(summary_data[start[j]:end[j]+1, 0]),
                 Rank     = np.mean(summary_data[start[j]:end[j]+1, 1]), 
                 n        = np.mean(summary_data[start[j]:end[j]+1, 2]), 
                 RAdeg    = np.mean(summary_data[start[j]:end[j]+1, 3]), 
                 DEdeg    = np.mean(summary_data[start[j]:end[j]+1, 4]), 
                 Rp       = np.mean(summary_data[start[j]:end[j]+1, 5]), 
                 Per      = np.mean(summary_data[start[j]:end[j]+1, 6]), 
                 S        = np.mean(summary_data[start[j]:end[j]+1, 7]), 
                 K        = np.mean(summary_data[start[j]:end[j]+1, 8]),
                 Rstar    = np.mean(summary_data[start[j]:end[j]+1, 9]), 
                 Teff     = np.mean(summary_data[start[j]:end[j]+1, 10]), 
                 Jmag     = np.mean(summary_data[start[j]:end[j]+1, 11]), 
                 a        = np.mean(summary_data[start[j]:end[j]+1, 12]),
                 Teq      = np.mean(summary_data[start[j]:end[j]+1, 13]),                   
                 gs       = np.mean(summary_data[start[j]:end[j]+1, 14]), 
                 Mmed     = np.mean(summary_data[start[j]:end[j]+1, 15]), 
                 Rhomed   = np.mean(summary_data[start[j]:end[j]+1, 18]), 
                 tdur     = np.mean(summary_data[start[j]:end[j]+1, 19]), 
                 nt10yr   = np.mean(summary_data[start[j]:end[j]+1, 20]),
                 nt_lo    = np.mean(summary_data[start[j]:end[j]+1, 21]),
                 nt_hi    = np.mean(summary_data[start[j]:end[j]+1, 22]),
                 tT_lo    = np.mean(summary_data[start[j]:end[j]+1, 23]),
                 tT_hi    = np.mean(summary_data[start[j]:end[j]+1, 24]),
                 tT_avg   = np.mean(summary_data[start[j]:end[j]+1, 25])))
    
            fout.write(fmt5.format(
                 Cat      = np.max(summary_data[start[j]:end[j]+1, 0]),
                 Rank     = np.max(summary_data[start[j]:end[j]+1, 1]), 
                 n        = np.max(summary_data[start[j]:end[j]+1, 2]), 
                 RAdeg    = np.max(summary_data[start[j]:end[j]+1, 3]), 
                 DEdeg    = np.max(summary_data[start[j]:end[j]+1, 4]), 
                 Rp       = np.max(summary_data[start[j]:end[j]+1, 5]), 
                 Per      = np.max(summary_data[start[j]:end[j]+1, 6]), 
                 S        = np.max(summary_data[start[j]:end[j]+1, 7]), 
                 K        = np.max(summary_data[start[j]:end[j]+1, 8]),
                 Rstar    = np.max(summary_data[start[j]:end[j]+1, 9]), 
                 Teff     = np.max(summary_data[start[j]:end[j]+1, 10]), 
                 Jmag     = np.max(summary_data[start[j]:end[j]+1, 11]), 
                 a        = np.max(summary_data[start[j]:end[j]+1, 12]),
                 Teq      = np.max(summary_data[start[j]:end[j]+1, 13]),                   
                 gs       = np.max(summary_data[start[j]:end[j]+1, 14]), 
                 Mmed     = np.max(summary_data[start[j]:end[j]+1, 15]), 
                 Rhomed   = np.max(summary_data[start[j]:end[j]+1, 18]), 
                 tdur     = np.max(summary_data[start[j]:end[j]+1, 19]), 
                 nt10yr   = np.max(summary_data[start[j]:end[j]+1, 20]),
                 nt_lo    = np.max(summary_data[start[j]:end[j]+1, 21]),
                 nt_hi    = np.max(summary_data[start[j]:end[j]+1, 22]),
                 tT_lo    = np.max(summary_data[start[j]:end[j]+1, 23]),
                 tT_hi    = np.max(summary_data[start[j]:end[j]+1, 24]),
                 tT_avg   = np.max(summary_data[start[j]:end[j]+1, 25])))
            
            break
            
             
# Compute stats for full target set    
fout.write('\n')    
fout.write('\n')
fout.write('\n')
fout.write('\n')
fout.write(head5)
fout.write(head1)
fout.write(head4)
fout.write(fmt6.format(
                     Cat       = np.min(summary_data[:, 0]),
                     Rank      = np.min(summary_data[:, 1]), 
                     n         = np.min(summary_data[:, 2]), 
                     RAdeg     = np.min(summary_data[:, 3]), 
                     DEdeg     = np.min(summary_data[:, 4]), 
                     Rp        = np.min(summary_data[:, 5]), 
                     Per       = np.min(summary_data[:, 6]), 
                     S         = np.min(summary_data[:, 7]), 
                     K         = np.min(summary_data[:, 8]),
                     Rstar     = np.min(summary_data[:, 9]), 
                     Teff      = np.min(summary_data[:, 10]), 
                     Jmag      = np.min(summary_data[:, 11]), 
                     a         = np.min(summary_data[:, 12]),
                     Teq       = np.min(summary_data[:, 13]),                   
                     gs        = np.min(summary_data[:, 14]), 
                     Mmed      = np.min(summary_data[:, 15]), 
                     Rhomed    = np.min(summary_data[:, 18]), 
                     tdur      = np.min(summary_data[:, 19]), 
                     nt10yr    = np.min(summary_data[:, 20]),
                     nt_lo     = np.min(summary_data[:, 21]),
                     nt_hi     = np.min(summary_data[:, 22]),
                     tT_lo     = np.min(summary_data[:, 23]),
                     tT_hi     = np.min(summary_data[:, 24]),
                     tT_avg    = np.min(summary_data[:, 25])))

fout.write(fmt7.format(
                     Cat       = np.median(summary_data[:, 0]),
                     Rank      = np.median(summary_data[:, 1]), 
                     n         = np.median(summary_data[:, 2]), 
                     RAdeg     = np.median(summary_data[:, 3]), 
                     DEdeg     = np.median(summary_data[:, 4]), 
                     Rp        = np.median(summary_data[:, 5]), 
                     Per       = np.median(summary_data[:, 6]), 
                     S         = np.median(summary_data[:, 7]), 
                     K         = np.median(summary_data[:, 8]),
                     Rstar     = np.median(summary_data[:, 9]), 
                     Teff      = np.median(summary_data[:, 10]), 
                     Jmag      = np.median(summary_data[:, 11]), 
                     a         = np.median(summary_data[:, 12]),
                     Teq       = np.median(summary_data[:, 13]),                   
                     gs        = np.median(summary_data[:, 14]), 
                     Mmed      = np.median(summary_data[:, 15]), 
                     Rhomed    = np.median(summary_data[:, 18]), 
                     tdur      = np.median(summary_data[:, 19]), 
                     nt10yr    = np.median(summary_data[:, 20]),
                     nt_lo     = np.median(summary_data[:, 21]),
                     nt_hi     = np.median(summary_data[:, 22]),
                     tT_lo     = np.median(summary_data[:, 23]),
                     tT_hi     = np.median(summary_data[:, 24]),
                     tT_avg    = np.median(summary_data[:, 25])))       
        
fout.write(fmt8.format(
                     Cat       = np.mean(summary_data[:, 0]),
                     Rank      = np.mean(summary_data[:, 1]), 
                     n         = np.mean(summary_data[:, 2]), 
                     RAdeg     = np.mean(summary_data[:, 3]), 
                     DEdeg     = np.mean(summary_data[:, 4]), 
                     Rp        = np.mean(summary_data[:, 5]), 
                     Per       = np.mean(summary_data[:, 6]), 
                     S         = np.mean(summary_data[:, 7]), 
                     K         = np.mean(summary_data[:, 8]),
                     Rstar     = np.mean(summary_data[:, 9]), 
                     Teff      = np.mean(summary_data[:, 10]), 
                     Jmag      = np.mean(summary_data[:, 11]), 
                     a         = np.mean(summary_data[:, 12]),
                     Teq       = np.mean(summary_data[:, 13]),                   
                     gs        = np.mean(summary_data[:, 14]), 
                     Mmed      = np.mean(summary_data[:, 15]), 
                     Rhomed    = np.mean(summary_data[:, 18]), 
                     tdur      = np.mean(summary_data[:, 19]), 
                     nt10yr    = np.mean(summary_data[:, 20]),
                     nt_lo     = np.mean(summary_data[:, 21]),
                     nt_hi     = np.mean(summary_data[:, 22]),
                     tT_lo     = np.mean(summary_data[:, 23]),
                     tT_hi     = np.mean(summary_data[:, 24]),
                     tT_avg    = np.mean(summary_data[:, 25])))        

fout.write(fmt9.format(
                     Cat       = np.max(summary_data[:, 0]),
                     Rank      = np.max(summary_data[:, 1]), 
                     n         = np.max(summary_data[:, 2]), 
                     RAdeg     = np.max(summary_data[:, 3]), 
                     DEdeg     = np.max(summary_data[:, 4]), 
                     Rp        = np.max(summary_data[:, 5]), 
                     Per       = np.max(summary_data[:, 6]), 
                     S         = np.max(summary_data[:, 7]), 
                     K         = np.max(summary_data[:, 8]),
                     Rstar     = np.max(summary_data[:, 9]), 
                     Teff      = np.max(summary_data[:, 10]), 
                     Jmag      = np.max(summary_data[:, 11]), 
                     a         = np.max(summary_data[:, 12]),
                     Teq       = np.max(summary_data[:, 13]),                   
                     gs        = np.max(summary_data[:, 14]), 
                     Mmed      = np.max(summary_data[:, 15]), 
                     Rhomed    = np.max(summary_data[:, 18]), 
                     tdur      = np.max(summary_data[:, 19]), 
                     nt10yr    = np.max(summary_data[:, 20]),
                     nt_lo     = np.max(summary_data[:, 21]),
                     nt_hi     = np.max(summary_data[:, 22]),
                     tT_lo     = np.max(summary_data[:, 23]),
                     tT_hi     = np.max(summary_data[:, 24]),
                     tT_avg    = np.max(summary_data[:, 25]))) 


fout.write(textwrap.dedent("""\


                           
Notes: 
    "995 and 9950.0"   tT_avg undetermined since hi_metal atm not detected 
    "996 and 9960.0"   Rp > 10:  Mmed estimates are unreliable
    "997 and 9970.0"   number of transits needed for detection exceeds"""
    + """ those observable in 10 yr mission
    "998 and 9980.0"   target Jmag is below (brighter than) the saturation"""
    + """ threshold of this instrument/mode
    "999 and 9990.0"   dBIC is < dBIC threshold for this target: """  
    + """no detection for this atm assumption    
    
 """))
    
fout.close()


print(' JET run complete!') 

raise SystemExit