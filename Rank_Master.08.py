

"""
@author: Charles Fortenbach

Code associated with paper submitted to AAS Journals (6/10/19)

Title: A FRAMEWORK FOR OPTIMIZING EXOPLANET TARGET SELECTION FOR
THE JAMES WEBB SPACE TELESCOPE

Subprogram:  Rank_Master.py

Version:  1.0

"""

import numpy as np
import textwrap
import datetime
import time
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
# See below

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
Nrow = Nrow_end - Nrow_start + 1
start = np.zeros(7, dtype=int)
end = np.zeros(7, dtype=int)



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
now = datetime.datetime.now()
timestr = time.strftime("%Y%m%d-%H%M%S")

fout = open(os.getcwd() + '/JET_Smry_tables/' +'JET_Smry_' 
            + timestr + '.txt', "w", newline='')
fout.write('JET run ID:  ' + str(now))
fout.write(textwrap.dedent("""\


 Target Summary Table:
                      
 Col  Label     Explanations
 --------------------------------------------------------------------------
 ---  ranking:
   0  Cat       Planet demographic category (1 - 7)
   1  Rank      Rank within Cat based on total observ. cycle for detection
  
 ---  survey data:
   2  Targ      Survey table index (starting with 1)
   3  RAdeg     Right Ascension in decimal degrees (J2000)
   4  DEdeg     Declination in decimal degrees (J2000)
   5  Rp        Planetary radius (mean) in Earth units (R/Rearth)
   6  Per       Period (days)
   7  S         Planetary insolation in Earth units (S/Searth)
   8  R*        Stellar radius (Solar radii)
   9  Teff      Stellar effective temperature (K)
  10  Jmag      Apparent J band magnitude (mag)
  
 ---  computed values:
  11  a         Semi-major axis of planet orbit (Solar radii)
  12  Teq       Planet equilibrium temperature (K)
  13  gs        Planet surface gravity (m/s^2)
  14  Mmed      Planet mass (median) in Earth units (M/Mearth)
  15  tdur      Transit duration (hrs)
  16  nt10yr    Number of transits observable in 10yr mission
  17  nt_lo     Number of transits needed for detection; lo-metal atm
  18  nt_hi     Number of transits needed for detection; hi-metal atm
  19  tT_lo     Total observ. cycle (in/out of transit); lo-metal atm (hrs)
  20  tT_hi     Total observ. cycle (in/out of transit); hi-metal atm (hrs)
  21  tT_avg    Total observ. cycle for 'average' atm detection (hrs)

 """))
    
fout.write(' JWST Instrument: ' + Instrument + '\n')  
fout.write(' Wavelength Range: ' + str(wavelimlo) + ' microns - ' 
           + str(wavelimhi) + ' microns' + '\n')  
fout.write(' Jmag Limit: ' + str(Jmag_sat) + '\n') 
fout.write(' Noise Floor: ' + str(int(nfloor)) + ' ppm' + '\n') 

if Cloud_lo == 0.0:
    fout.write(' Eq. of State (lo_metal atm): ' + str(EOS_lo) 
    + ',  with no clouds' + '\n')    
else:
    fout.write(' Eq. of State (lo_metal atm): ' + str(EOS_lo) 
    + ',  with cloud deck at ' + str(Cloud_lo/100) + ' mbar' + '\n')
    
if Cloud_hi == 0.0:
    fout.write(' Eq. of State (hi_metal atm): ' + str(EOS_hi) 
    + ',  with no clouds' + '\n')    
else:
    fout.write(' Eq. of State (hi_metal atm): ' + str(EOS_hi) 
    + ',  with cloud deck at ' + str(Cloud_hi/100) + ' mbar' + '\n') 
    
fout.write(' Detection Threshold (dBIC): ' + str(dBIC_thresh) + '\n')  

    

# Define headers and formats for output
head1 = ("      |------------------------------survey data---------------------------|"
         "-------------------------------------computed values---------------------------|\n")
head2 = (" Rank   Targ    RAdeg    DEdeg     Rp    Per       S     R*   Teff   "
         "Jmag      a    Teq     gs    Mmed   tdur  nt10yr  nt_lo  nt_hi   "
         "tT_lo   tT_hi  tT_avg\n")
fmt1 = ("{Rank:5.0f} {n:6.0f} {RAdeg:8.3f} {DEdeg:8.3f} {Rp:6.2f} {Per:6.2f}"
        " {S:7.1f} {Rstar:6.2f} {Teff:6.0f} {Jmag:6.1f} {a:6.1f} {Teq:6.0f}"
        " {gs:6.1f} {Mmed:7.2f} {tdur:6.1f} {nt10yr:7.0f} {nt_lo:6.0f}"
        " {nt_hi:6.0f} {tT_lo:7.1f} {tT_hi:7.1f} {tT_avg:7.1f}\n")

head3 = ("        Targ    RAdeg    DEdeg     Rp    Per       S     R*   Teff   "
         "Jmag      a    Teq     gs    Mmed   tdur  nt10yr\n")
fmt2 = (" Min:                          {Rp:6.2f} {Per:6.2f} {S:7.1f}"
        " {Rstar:6.2f} {Teff:6.0f} {Jmag:6.1f} {a:6.1f} {Teq:6.0f} {gs:6.1f}"
        " {Mmed:7.2f} {tdur:6.1f} {nt10yr:7.0f}\n")
fmt3 = (" Med:                          {Rp:6.2f} {Per:6.2f} {S:7.1f}"
        " {Rstar:6.2f} {Teff:6.0f} {Jmag:6.1f} {a:6.1f} {Teq:6.0f} {gs:6.1f}"
        " {Mmed:7.2f} {tdur:6.1f} {nt10yr:7.0f}\n")
fmt4 = ("Mean:                          {Rp:6.2f} {Per:6.2f} {S:7.1f}"
        " {Rstar:6.2f} {Teff:6.0f} {Jmag:6.1f} {a:6.1f} {Teq:6.0f} {gs:6.1f}"
        " {Mmed:7.2f} {tdur:6.1f} {nt10yr:7.0f}\n")
fmt5 = (" Max:                          {Rp:6.2f} {Per:6.2f} {S:7.1f}"
        " {Rstar:6.2f} {Teff:6.0f} {Jmag:6.1f} {a:6.1f} {Teq:6.0f} {gs:6.1f}"
        " {Mmed:7.2f} {tdur:6.1f} {nt10yr:7.0f}\n")

head5 = "Full Target Set:\n" 
head4 = ("        Targ    RAdeg    DEdeg     Rp    Per       S     R*   Teff   "
         "Jmag      a    Teq     gs    Mmed   tdur  nt10yr\n")
fmt6 = (" Min: {n:6.0f} {RAdeg:8.3f} {DEdeg:8.3f} {Rp:6.2f} {Per:6.2f}"
        " {S:7.1f} {Rstar:6.2f} {Teff:6.0f} {Jmag:6.1f} {a:6.1f} {Teq:6.0f}"
        " {gs:6.1f} {Mmed:7.2f} {tdur:6.1f} {nt10yr:7.0f}\n")
fmt7 = (" Med:        {RAdeg:8.3f} {DEdeg:8.3f} {Rp:6.2f} {Per:6.2f} {S:7.1f}"
        " {Rstar:6.2f} {Teff:6.0f} {Jmag:6.1f} {a:6.1f} {Teq:6.0f} {gs:6.1f}"
        " {Mmed:7.2f} {tdur:6.1f} {nt10yr:7.0f}\n")
fmt8 = ("Mean:        {RAdeg:8.3f} {DEdeg:8.3f} {Rp:6.2f} {Per:6.2f} {S:7.1f}"
        " {Rstar:6.2f} {Teff:6.0f} {Jmag:6.1f} {a:6.1f} {Teq:6.0f} {gs:6.1f}"
        " {Mmed:7.2f} {tdur:6.1f} {nt10yr:7.0f}\n")
fmt9 = (" Max: {n:6.0f} {RAdeg:8.3f} {DEdeg:8.3f} {Rp:6.2f} {Per:6.2f} {S:7.1f}"
        " {Rstar:6.2f} {Teff:6.0f} {Jmag:6.1f} {a:6.1f} {Teq:6.0f} {gs:6.1f}"
        " {Mmed:7.2f} {tdur:6.1f} {nt10yr:7.0f}\n")


# Set start and end pts for categories
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

j = np.max(np.int(np.int(summary_data[i, 0])))-1
for i in range(0, Nrow_end - Nrow_start + 1): 
    if i == (Nrow_end - Nrow_start): 
        end[j] = i


# Set flag for targets with tT_avg undetermined since hi_metal atm not detected
for i in range(0, Nrow_end - Nrow_start + 1):
    if summary_data[i, 25] > 4000 and summary_data[i, 25] < 6000:
        summary_data[i, 25] = 9950    


# Write summary table
for j in range(0, 7):
    fout.write('\n')
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
fout.write(head3)
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
    "995 and 9950"   tT_avg undetermined since hi_metal atm not detected 
    "996 and 9960"   Rp > 10:  Mmed estimates are unreliable
    "997 and 9970"   number of transits needed for detection exceeds those"""
    + """ observable in 10 yr mission
    "998 and 9980"   target Jmag is below (brighter than) the saturation"""
    + """ threshold of this instrument/mode
    "999 and 9990"   dBIC is < dBIC threshold for this target:  no detection"""
    + """ for this atm assumption    
    
 """))
    
fout.close()


print(' JET run complete!') 

raise SystemExit