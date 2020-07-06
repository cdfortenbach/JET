# *JWST* Exoplanet Targeting program (JET)
Computer tool to rank lists of exoplanet targets for atmospheric characterization by *JWST*.

## Description
The *James Webb Space Telescope* (*JWST*) will devote significant observing time to the study of exoplanets. It will not be serviceable as was the Hubble Space Telescope, and therefore the spacecraft/instruments will have a relatively limited life.  It is important to get as much science as possible out of this limited observing time.  We provide a computer tool that can be used to optimize lists of exoplanet targets for atmospheric characterization. Our tools take catalogs of planet detections, either simulated, or actual; categorize the targets by radius and equilibrium temperature; estimate planet masses; generate model spectra and simulated instrument spectra; perform a statistical analysis to determine if the instrument spectra can confirm an atmospheric detection; and finally, rank the targets within each category by observation time required for detection.

## Getting Started

### Hardware Configuration
* The JET code was developed and tested on a modest custom-built desktop computer (see Table 1).  In addition, all production runs to date have been executed on this machine.  

 Table 1. Hardware for Development and Testing

| Component         | Details                |
|:------------------|:-----------------------|
| Processor         | Intel Core i7-2600     |
| Clock             | 3.40 GHz               |
| CPU's             | 4                      |
| Threads/CPU       | 2                      |
| Motherboard       | ASUSTeK P8Z68-V        |
| RAM               | 16 GB                  |
| Internal Storage  | 256 GB SSD             |
| External Storage  | 3 TB                   |
| Cooling           | Air                    |


* In hindsight this computer system was fine for development, but for test runs with the large target survey data-set the processing speed was inadequate.  The overall run time was ~ 3.4 minutes per target (in a single processor mode).  With just under 2000 targets in the baseline survey this meant a total run time of ~ 110 hours.  We ultimately decided to split the target list and run two instantiations of the code in parallel.  This reduced the overall run time to ~ 55 hours.

* An effort to make the code more efficient would be helpful. In particular, consideration of processing multiple targets with a parallel processing architecture would seem to be a worthwhile development project; however, this idea has not been further addressed here.

* The full JET code has a memory storage footprint of approximately 3.7 GB.  A production run of a 2000 target catalog will generate about 1 GB of output (spectra, and various tables).  So, you need to plan on a total memory storage requirement of at least 4.7 GB for a large survey run.  In execution you will need about 4 GB of free RAM on top of any other system requirements.

### Software Configuration/Installation
* The desktop computer described in the previous section runs the Windows 7 Ultimate (SP1) operating system.  We chose to develop the JET suite of programs to run in a Linux environment.  We used Oracle VirtualBox software (5.1.26) to run a Linux guest OS (specifically Ubuntu 16.04) under the Windows 7 host.  Of course a dedicated Linux machine should work just fine.  The installation and execution instructions provided here assume that the user will be using a Linux OS.

* The VirtualBox software allows files to be shared between the two OS's, which is very useful, but otherwise it keeps the two systems separate.  In order to implement the file sharing feature it was necessary to edit the system BIOS to enable the CPU's virtualization features.  Unfortunately, the details of the editing are different with different motherboards and CPU hardware so it is hard to give specifics here.  The key is to sort through the BIOS menus at boot up to find the CPU virtualization feature toggle and enable it.  Of course, if the JET suite is installed on a Linux machine in the first place there is no need to worry about virtualization software.

* The JET code was written in Python 3.6.1/2.  It is highly recommended that you use Python 3.6.2 to run JET.

* JET has a number of Python dependencies (see Table 2); however, many of them will be satisfied by an up-to-date scientific Python installation.  Unless you actively maintain your own scientific Python distribution, we recommend installing the latest Anaconda Linux distribution (currently, as of 7-5-2019, for Python 3.7) and running the following command to downgrade to Python 3.6 in the root environment:

 ```conda install python=3.6.2```

* To avoid conflicts, it is recommended that you first remove/uninstall any existing packages shown in Table 2 with the incorrect version/build.  You can see the packages available in your Anaconda distribution with the command:  ```conda list```.

* To remove packages that do not have the < pip > designation, you can use the command:

 ```conda remove package_name```

* For packages with the < pip > designation, you can use the command:

 ```pip uninstall package_name```

* **Warning!**:  rolling back certain elements of your Python package distribution, may break other Python programs that are dependent on your existing package distribution.  It may make sense to use virtualenv (https://pypi.org/project/virtualenv/), a tool to create isolated Python environments.

* For those packages that are not included in the basic conda distribution, you will need to enable certain additional distribution channels with:

 ```conda config --add channels conda-forge```

 and

 ```conda config --add channels http://ssb.stsci.edu/astroconda```

 This only needs to be done once, not for each package.

* It is recommended that you install the version/build of the packages specified in Table 2. There are various complex dependencies and this build has been verified.  

* Many of the required packages (unless they have < pip > in the build column of Table 2) can then be installed from the Linux command line with:

 ```conda install package_name=version=build_string```

 For example:

 ```conda install joblib=0.11=py36_0```

* For those packages with a non-conda channel, you should use the command:

 ```conda install -c channel package_name=version=build_string```

 For example:

 ```conda install -c conda-forge astropy=3.2.3```

 or another example:

 ```conda install -c http://ssb.stsci.edu/astroconda photutils=0.4.1```

* For the remaining packages, you will need to use the pip package manager.  For example, to install numpy use the command:

 ```pip install numpy=1.14.0```

* When this process has been completed for all of the required packages, again use the Linux command ```conda list``` to verify that all packages and modules shown in Table 2 are installed:

 Table 2. Python Packages/Modules for JET Installation

| Name                  | Type              | Version/Build                          | Channel                             |
|:----------------------|:------------------|:---------------------------------------|:------------------------------------|
| astropy               | Python pkg.       | 3.2.3                                  |                                     |
| bokeh                 | Python pkg.       | 0.12.6 / py36_0                        | conda-forge                         |
| joblib                | Python pkg.       | 0.11 / py36_0                          | conda-forge                         |
| matplotlib            | Python pkg.       | 2.2.2                                  |                                     |
| numpy                 | Python pkg.       | 1.14.0 / < pip >                       | pypi                                |
| pandeia.engine        | Python pkg.       | 1.4 / < pip >                          | pypi                                |
| pandexo.engine        | Python pkg.       | 1.3 / < pip >                          | pypi                                |
| photutils             | Python pkg.       | 0.4.1                                  | http://ssb.stsci.edu/astroconda     |
| pyephem               | Python pkg.       | 3.7.6.0 / < pip >                      | pypi                                |
| pyfftw                | Python pkg.       | 0.10.4                                 | conda-forge                         |
| pysynphot             | Python pkg.       | 0.9.8.8 / py36_1                       | http://ssb.stsci.edu/astroconda     |
| scipy                 | Python pkg.       | 1.0.0 / py36_blas_openblas_201         | conda-forge                         |
| spectres              | Python module     | 2.0.0 / < pip >                        | pypi                                |
| sphinx                | Python pkg.       | 1.5.6                                  | conda-forge                         |
| synphot               | Python pkg.       | 0.1.3                                  | http://ssb.stsci.edu/astroconda     |

### Installing the JET Code

1. Now, choose/create a working directory for the JET code within your Linux system file structure.

2. Open a terminal in your Linux Download directory.  Then, download the JET repository tar file from GitHub with the command:

 ```curl -L https://github.com/cdfortenbach/JET/tarball/master > master```

   This is a big file, so it may take a few minutes.  You should eventually see a folder/directory called **master** appear in the Download directory.

3. Now, unpack the **master** tar file with the command:

 ```tar -xvzf master```

   you should now see a new folder/directory called **cdfortenbach-JET-*commitID#***.   

4. Open this folder/directory and move the contents to your working directory.  

5. The downloaded tar file (**master**), and the now empty folder/directory **cdfortenbach-JET-*commitID#*** can be deleted.  

6. You should verify that you have all files and folders listed in Table 3.

 Table 3. Files and Folders in JET Working Directory

| Files/Folders               | Details                                                                                  |
|:----------------------------|:-----------------------------------------------------------------------------------------|
| dBIC                        | empty folder/directory to collect dBIC plots                                             |
| EOS                         | folder of equation of state data files for Exo-Transmit                                  |
| JET_Smry_tables             | empty folder/directory to collect Smry tables                                            |
| Opac                        | folder of opacity data files for Exo-Transmit                                            |
| Pdxo_Output_archive         | empty folder/directory to collect Pdxo_Output.txt files                                  |                 
| Spectra                     | empty folder/directory to collect spectral plots and tables                              |
| T_P                         | folder of temperature-pressure profile files for Exo-Transmit                            |
| Target_Surveys              | folder of various example target surveys in format of "Sullivan" survey .mrt             |
| Append_util.py              | Python utility associated with multiple instantiations of JET                            |
| Exo_Transmit                | atmospheric modelling code (C program, compiled-executable)                              |   
| ExoT_Master.59.py           | Python program (element of JET)                                                          |
| ExoT_USR_MANUAL.txt         | USR_MANUAL for Exo-Transmit                                                              |
| JET_Input.txt               | file for primary input to JET                                                            |
| JET_Master.57.sh            | Linux shell script controlling the JET code                                              |
| LICENSE                     | Open-source license GNU GPLv3                                                            |
| otherInput.in               | other input information for Exo_Transmit                                                 |
| Pdxo_Input.txt              | internal data transfer file                                                              |
| Pdxo_Master.69.py           | Python program (element of JET)                                                          |
| Pdxo_Output.txt             | internal data transfer file                                                              |
| Rank_Master.11.py           | Python program (element of JET)                                                          |
| README.md                   | README for JET in markdown format                                                        |
| selectChem.in               | chemistry file selector system for Exo-Transmit                                          |
| target_survey.txt           | .txt file in format of the "Sullivan" survey .mrt                                        |
| userInput.in                | internal transfer file assoc. with Exo-Transmit                                          |


7. Next, you need to download certain data-files associated with PandExo/Pandeia.  The first is a set of files for Pandeia itself, which can be downloaded from: https://stsci.app.box.com/v/pandeia-refdata-v1p4.  Once the download is complete, move the tar file from your download folder to the JET working directory.  This is in gzip form and can be unpacked in the working directory. Use the command:

 ```$ tar -xvzf pandeia_data-1.4.tar.gz```  

  to unpack the tar file. Once the files are unpacked you should be left with a folder/directory named **pandeia_data-1.4**.  The tar file can then be deleted.

8. Next, you will need to download another set of data-files associated with PandExo.  The second set of data-files are for pysynphot. Pandeia uses pysynphot internally for creating reference spectra. The pysynphot reference files may be downloaded from: https://archive.stsci.edu/pub/hst/pysynphot/.  We want to download **synphot5.tar.gz**.  Move the tar file from your download folder to the JET working directory.  This is also in gzip form and can be unpacked in the working directory with the command:

 ```$ tar -xvzf synphot5.tar.gz```

  Once the files are unpacked you should be left with a folder/directory named **grp** and the tar file.  The tar file can be deleted.

 Your working directory should now have all of the necessary elements of JET, including all files and folders associated with Exo-Transmit and PandExo. Of course it is possible to download the latest version of Exo-Transmit (source code) from GitHub (see https://github.com/elizakempton/Exo_Transmit), but for simplicity we have chosen to include a fully compiled, and tested version of Exo-Transmit as part of the JET download.  There are also some minor changes to a couple of the baseline Exo-Transmit data-files that we have made for JET to work properly.

  In some cases, high-metallicity atmospheres in particular, can be optically thick all the way to the very top of the atmosphere, at certain wavelengths.  This will cause a hard cutoff of the spectrum at a specific transit depth.  To extend the atmosphere to higher levels (i.e., lower pressures), more lines have been added to the original release T-P files (we went from 333 to 533 lines), extending the decaying exponential to lower pressure.  The number of optical depth points given in the file: **otherInput.in**, has also been be modified to give the correct number of lines for the new T-P files.  You should not need to edit these files further.

  There are several input files associated with Exo-Transmit that have not been described in detail (e.g., **selectChem.in**, **otherInput.in**, and **userInput.in**).  For our purposes these should not have to be disturbed by the user; however, it is possible to make changes to these if necessary.  Further details can be found in the Exo-Transmit user manual included in this repository, or at https://github.com/elizakempton/Exo_Transmit.

9. It is recommended that the user edit certain print statements in the **PandExo.engine** (installed as a Python package in the Anaconda directory) that are unnecessary and time consuming in our long, many-target runs. Specifically, we recommend that you "comment out" (add leading # character to the line) the following lines in **justdoit.py** and **jwst.py**:

    anconda3/lib/python3.6/site-packages/pandexo/engine/justdoit.py - line 267
    	should be:   ```#print("Running Single Case for: " + inst[0])```

    anconda3/lib/python3.6/site-packages/pandexo/engine/jwst.py - line 153
    	should be:   ```#print("Optimization Reqested: Computing Duty Cycle")```

    anconda3/lib/python3.6/site-packages/pandexo/engine/jwst.py - line 156
    	should be:   ```#print("Finished Duty Cycle Calc")```

    anconda3/lib/python3.6/site-packages/pandexo/engine/jwst.py - line 162
    	should be:   ```#print("Starting Out of Transit Simulation")```

    anconda3/lib/python3.6/site-packages/pandexo/engine/jwst.py - line 169
    	should be:   ```#print("End out of Transit")```

    anconda3/lib/python3.6/site-packages/pandexo/engine/jwst.py - line 176
    	should be:   ```#print("Starting In Transit Simulation")```

    anconda3/lib/python3.6/site-packages/pandexo/engine/jwst.py - line 178
    	should be:   ```#print("End In Transit")```

10. That’s it!  You have completed the JET installation.


### Executing the program

The following step-by-step process should guide you through making a run with the JET code:

1.  First, using a GUI (e.g., Ubuntu, etc.) or with a Linux terminal, navigate to the JET working directory.

2.  If you want to start fresh, then the folders **dBIC**, **Spectra**, **Pdxo_Output_archive**, and **JET_Smry_tables** should be cleaned out for a new run.  Any files in these folders/directories can be saved elsewhere as appropriate. The dummy **README** file in each of these folder/directories can be left in place.

3.  Using a .txt editor (or development environment), open the **JET_Input.txt** file and make the appropriate changes to the input parameters.  The following shows an example set-up for a run with NIRSpec G395M.  Only the boldface parameters should be changed.  Only make changes to User Input values in rows 7 - 46, and cols 60 - 79!!   Do not change the syntax for any item.  This is a comma delimited file.  Commas should only be included where indicated.

   JET_Input.txt file:

| Parameter Description                                    | User Input                |
|:---------------------------------------------------------|:--------------------------|
| Starting target from catalog for this run:               |, **5**|   
| Ending target from catalog:                              |, **7**|
|                                                          |                           |
| JWST instrument:                                         |, **NIRSpec G395M**|
| Wavelength short limit (microns):                        |, **2.87**|
| Wavelength long limit (microns):                         |, **5.18**|
| Jmag limit (Teff = 10000K):                              |, **6.2**|
| Jmag limit (Teff = 5000K):                               |, **6.8**|
| Jmag limit (Teff = 2500K):                               |, **7.4**|
| Detector linear response limit (% FW):                   |, **80**|
| Noise floor (nfloor) (ppm):                              |, **25**|
| R value of sim (Res):                                    |, **100**|
| Number of dBIC samples for each ntr grid pt.:            |, **2000**|
| Detection threshold (dBIC):                              |, **10**|
| Free model spectrum BIC parameters:                      |, **5**|
|                                                          |                           |
| Eq. of State (lo metal atm):                             |, **eos_5Xsolar_cond**|
| Cloud_lo (Pa):                                           |, **0**|
| Eq. of State (hi metal atm):                             |, **eos_1000Xsolar_cond**|
| Cloud_hi (Pa):                                           |, **10000**|
|                                                          |                           |
| Out of transit factor (% tdur):                          |, **100**|
| \+ Out of transit "timing tax" (sec):                    |, **6300**|
| Slew duration avg. (sec):                                |, **1800**|
| SAMs: small angle maneuvers (sec):                       |, **138**|
| GS Acq: guide star acquisition(s)(sec):                  |, **312**|
| Targ Acq: target acquisition if any (sec):               |, **62**|
| Exposure Ovhd: factor 1:                                 |, **0.0536**|
| Exposure Ovhd: factor 2 (sec):                           |, **0**|
| Mech: mechanism movements (sec):                         |, **220**|
| OSS: Onboard Script System compilation (sec):            |, **58**|
| MSA: NIRSpec MSA configuration (sec):                    |, **0**|
| IRS2: NIRSpec IRS2 Detector Mode setup (sec):            |, **217**|
| Visit Ovhd: visit cleanup activities (sec):              |, **102**|
| Obs Ovhd factor (%):                                     |, **16**|
| DS Ovhd (sec):                                           |, **0**|        
|                                                          |                           |
| RunExoT (Y/N):                                           |, **Y**|
| RunPdxo (Y/N):                                           |, **Y**|
| RunRank (Y/N):                                           |, **Y**|


  * The starting and ending rows to analyze from the target_survey.txt file need to be specified.

  * The instrument/mode (e.g., NIRSpec G395M) should be selected from the listing (comments further down in the .txt file).  Four of the modes have been tested.  

  * The wavelength ranges for the given instrument/mode need to be set.  The ranges for the various instrument/modes are shown further down in the .txt file comments.  

  * The various instrument/modes have brightness limits that must be set. The "full-well" limits are shown further down in the .txt file comments for the particular instrument/mode studied.  The detector linearity limit as a % of "full-well" must also be set.

  * The multi-transit residual noise-floor for the given instrument/mode needs to be set.  Suggested values are shown further down in the .txt file comments.

  * The spectral resolving power (R value = Res) should be set.  The model spectra are initially generated with a resolving power of ~ 1000.  For the statistical analysis and for plotting purposes the model spectra and simulated spectra are binned down to the same resolving power.  The factor, Res, sets this value.  This is generally well below the native resolving power of the instrument/modes studied.

  * The number of random noise samples for each of the dBIC vs number-of-transit grid points must be set.  We believe that an appropriate compromise between computation time and result stability is 2000 samples.

  * We are generally using a dBIC_threshold of 10 in our BIC model selection procedure.  This value should give us a "very strong" detection.

  * We have set the number of free model spectrum BIC parameters to five for our Baseline.

  * The equation of state files (EOS\_lo and EOS\_hi), and cloud top pressure levels (Cloud\_lo and Cloud\_hi) also need to be set.  Further details on the EOS file names and cloud top pressure can be found in the Exo-Transmit user documentation (e.g., **ExoT\_USR_MANUAL.txt** included in the JET repository).

  * The various timing overheads are described in the associated paper.  The comments provide suggested values for these parameters.  In general these will only need to be changed for a change of instrument.

  * Sometimes it is useful to bypass certain parts of the code.  For example, if the model atmospheres have already been generated for a batch of targets.  It may makes sense to bypass the ExoT_Master section of the code.  The parameters RunExoT, RunPdxo, and RunRank can be toggled Y or N, to indicate whether or not to bypass that code element.  For a full run we will use all three of these parameters set to Y.


4. Next, make sure that the **target\_survey.txt** file in the working directory has the appropriate number of rows for your run.  Various example **target\_survey.txt** files are included in the **Target\_Survey** folder/directory in the repository.  
  * For our example run we will use the full Sullivan survey/catalog **target\_survey\_(Sullivan_survey).txt**.  The **target\_survey.txt** file in the repository is the full Sullivan survey/catalog so for our initial example run you should not need to make changes.

5. Next, execute the program.  From the Linux command line (from the working directory) enter the following:
 ```
 ./JET_Master.57.sh
 ```

   * After a few seconds, you should see . . .

 ```dos
 !Warning: The starting target value (Nrow_start) is not in sequence.

  If you proceed you may overwrite the previous entries.

  Proceed with ExoT? (Y/N):
 ```
6. You should enter **Y** (for the example run we are assuming that you are starting from scratch).

   * After a few seconds, you should see . . .
 ```dos
 Running ExoT_Master:

 Reading in survey data . . .

 Computing model transmission spectra using Exo-Transmit:

      Generating model spectra for target:  5

      Generating model spectra for target:  6

      Generating model spectra for target:  7

 Elapsed time for ExoT_Master (sec):  130.65741515159607

 !Warning: The starting target value (Nrow_start) is not in sequence.

 If you proceed you may overwrite the previous entries.

 Proceed with Pdxo? (Y/N):
 ```

7. You should enter **Y** (for our example run we are assuming that we are starting from scratch).

   * After a short time, you should see . . .
```dos
 Running Pdxo_Master:

 Transferring planetary system data and model spectra to Pandexo . . .

      Generating simulated instrument spectra for target:  5

         Comparing simulation to model spectrum and flat line for lo_metal atm:

               min # of transits for mean dBIC (less one sigma) > 10 :  1


         Comparing simulation to model spectrum and flat line for hi_metal atm:

               mean dBIC (less one sigma) < 10 for 50 transits ---> no detection


      Generating simulated instrument spectra for target:  6

         Comparing simulation to model spectrum and flat line for lo_metal atm:

               min # of transits for mean dBIC (less one sigma) > 10 :  1


         Comparing simulation to model spectrum and flat line for hi_metal atm:

               min # of transits for mean dBIC (less one sigma) > 10 :  47


      Generating simulated instrument spectra for target:  7

         Comparing simulation to model spectrum and flat line for lo_metal atm:

               min # of transits for mean dBIC (less one sigma) > 10 :  1


         Comparing simulation to model spectrum and flat line for hi_metal atm:

               min # of transits for mean dBIC (less one sigma) > 10 :  11


 Elapsed time for Pdxo_Master (sec):  170.68561053276062


 Running Rank_Master:

 Sorting and ranking targets . . .

 JET run complete!
```
8.  Now we can examine the products of the run.  Two summary tables should have been deposited into the **JET\_Smry_tables** folder/directory.  They are time and date stamped, and named something like: **JET_Smry_20190922-093016_un-ranked**, and **JET_Smry_20190922-093016.txt**.

 We can also examine the spectral data tables, and plots that have been deposited into the **Spectra** folder/directory.   For our test run you should see the files listed in Table 4.  

 Table 4. Files deposited in **Spectra** folder associated with test run:

| File name                       | Description                                                          |
|:--------------------------------|:---------------------------------------------------------------------|
| Trans_Spec_ExoT_5_hi_metal.txt  | full res model spectrum (data only) for hi-metal atm, targ 5         |
| Trans_Spec_ExoT_5_lo_metal.txt  | full res model spectrum (data only) for lo-metal atm, targ 5         |
| Trans_Spec_ExoT_6_hi_metal.txt  | full res model spectrum (data only) for hi-metal atm, targ 6         |
| Trans_Spec_ExoT_6_lo_metal.txt  | full res model spectrum (data only) for lo-metal atm, targ 6         |
| Trans_Spec_ExoT_7_hi_metal.txt  | full res model spectrum (data only) for hi-metal atm, targ 7         |
| Trans_Spec_ExoT_7_lo_metal.txt  | full res model spectrum (data only) for lo-metal atm, targ 7         |
| Trans_Spec_Pdxo_5_hi_metal.svg  | plot of re-binned model/instrument spectrum for hi-metal atm, targ 5 |
| Trans_Spec_Pdxo_5_lo_metal.svg  | plot of re-binned model/instrument spectrum for lo-metal atm, targ 5 |
| Trans_Spec_Pdxo_6_hi_metal.svg  | plot of re-binned model/instrument spectrum for hi-metal atm, targ 6 |
| Trans_Spec_Pdxo_6_lo_metal.svg  | plot of re-binned model/instrument spectrum for lo-metal atm, targ 6 |
| Trans_Spec_Pdxo_7_hi_metal.svg  | plot of re-binned model/instrument spectrum for hi-metal atm, targ 7 |
| Trans_Spec_Pdxo_7_lo_metal.svg  | plot of re-binned model/instrument spectrum for lo-metal atm, targ 7 |
| Trans_Spec_Xfer_5_hi_metal.txt  | xfer data of re-binned model spectrum for hi-metal atm, targ 5       |
| Trans_Spec_Xfer_5_lo_metal.txt  | xfer data of re-binned model spectrum for lo-metal atm, targ 5       |
| Trans_Spec_Xfer_6_hi_metal.txt  | xfer data of re-binned model spectrum for hi-metal atm, targ 6       |
| Trans_Spec_Xfer_6_lo_metal.txt  | xfer data of re-binned model spectrum for lo-metal atm, targ 6       |
| Trans_Spec_Xfer_7_hi_metal.txt  | xfer data of re-binned model spectrum for hi-metal atm, targ 7       |
| Trans_Spec_Xfer_7_lo_metal.txt  | xfer data of re-binned model spectrum for lo-metal atm, targ 7       |

 The dBIC-vs-number-of-transits plots have been deposited into the **dBIC** folder/directory.  For our test run you should see only two files (**dBIC_6_hi_metal.svg** and **dBIC_7_hi_metal.svg** ).  Some of the plots may appear to be missing.  They are not really missing, but rather they were not plotted due to a non-detection condition (e.g., saturation, etc.).  

 A copy of the **Pdxo_Output.txt** file (with target limits noted **Pdxo_Output_(5 - 7).txt**) has been deposited into the **Pdxo_Output_archive** folder.  This is a copy of the latest appended version of the **Pdxo_Output.txt** file.  This is essentially a backup if anything happens to the raw **Pdxo_Output.txt** file (e.g. an inadvertent overwrite, etc.).   If you are doing a series of runs where the targets are sequential then this file is being appended for each run.  It is over-written if the targets are out of sequence.

 You have completed the example run.  Congratulations!

### Additional Tools

As we've mentioned, the overall JET run time is ~ 3.4 minutes per target (in a single processor mode).  If you are interested in running a large target survey/catalog like "Sullivan" through JET, it is possible to reduce the clock time required by duplicating the JET working directory and splitting the target list.  The process proceeds normally as one instantiation of JET runs the first set of the split list (say from target 1 - 950), and another instantiation of JET runs the other set of the split list (from 951 - 1984).

After the runs are completed, in order to tie the two lists together and provide a complete ranked summary table, you can use the simple tool **Append_util.py** that will append the two **Pdxo_Output.txt** result files.  This tool can be launched from the low numbered list working directory with the command:
 ```python Append_util.py```  

The tool will ask for the file structure path to the high numbered list working directory.  Then it will perform the append operation and create a combined **Pdxo_Output.txt** file.  It will also deposit the combined file **Pdxo_Output (1-1984).txt** into the **Pdxo_Output_archive** folder/directory.

Recent update (as of 04-12-2020): The latest overall JET run time using an Intel i9-9900 CPU, with 6 of 16 cores assigned to a Virtual machine, is ~ 1.5 minutes per target for NIRSpec.  The run time for NIRISS is considerably longer, ~ 7.3 minutes per target.

## Author
Charles Fortenbach

## Citations
When publishing results based on usage of JET please cite:

for JET:
Fortenbach, C. D., & Dressing, C. D., 2020, arXiv:2002.01495 (https://ui.adsabs.harvard.edu/abs/2020arXiv200201495F)

and the following dependencies:

Exo-Transmit:
Kempton, E. M.-R., Lupu, R. E., Owusu-Asare, A., Slough, P., & Cale, B., 2016, arXiv:1611.03871 (http://adsabs.harvard.edu/abs/2016arXiv161103871K)

Freedman, R. S., Marley, M. S., & Lodders, K., 2008, ApJS, 174, 504-513

Freedman, R. S., Lustig-Yaeger, J., Fortney, J. J., et al., 2014, ApJS, 214, 25

Lupu, R. E., Zahnle, K., Marley, M. S., et al., 2014, ApJ, 784, 27

PandExo:
Batalha, N. E., Mandell, A., Pontoppidan, K., et al. 2017, Publications of the Astronomical Society of the Pacific, 129, 064501

astropy:
The Astropy Collaboration, Price-Whelan, A. M., Sip}ocz, B. M., et al. 2018, ArXiv e-prints, arXiv:1801.02634

matplotlib:
Hunter, J. D. 2007, Computing In Science & Engineering, 9, 90

numpy:
Oliphant, T. E. 2015, Guide to NumPy, 2nd edn. (USA: CreateSpace Independent Publishing Platform)

scipy:
Jones, E., Oliphant, T., Peterson, P., et al. 2001, SciPy: Open source scientific tools for Python

SpectRes:
Carnall, A. C. 2017, arXiv:1705.05165, 1705.05165v1

## Version History
* 0.1.1 (2018-12-10)
    Initial public release

* 0.1.2 (2019-07-15)
    Upgrade of pandexo.engine to v1.3, and pandeia.engine to v1.4.
    Update of Pandeia data-files to pandeia_data-1.4.
    Made path values indirect to enable easier portability.
    Modified path to grp directory in shell script for simpler installation.

* 1.1.0 (2019-10-01)
    New routine to handle instrument brightness limits.
    Increased number of model parameters from 4 to 5 for BIC calculation.
    Provided more sophisticated approach to estimating observing time.
    Various minor corrections and updates.

* 1.2.0 (2019-12-04)
    Made the number of free model spectrum BIC parameters a user input.

* 1.3.0 (2020-04-12)
    Updated the calculation of observation time elements, t_dwell and Science_time in Pdxo_Master based on latest STScI guidance.
    Updated JET_Master to now call updated Pdxo_Master.
    Updated package dependency table (Table 2), and File/Folder list (Table 3) in README.
    Updated JET_Input.txt time observation elements for active table, and commented tables for both NIRSpec and NIRISS based on latest STScI guidance.  

## License
&copy; 2019 Charles Fortenbach

This project is licensed under the GNU GPLv3 License - see the LICENSE.md file for details.

## Acknowledgments
The code discussed here was developed as part of C. Fortenbach's Master's thesis project at San Francisco State University. Prof. Courtney Dressing (UC Berkeley) was the thesis advisor on this effort and provided guidance and inspiration throughout.  C.F. would like to acknowledge the support of Prof. Mohsen Janatpour (College of San Mateo), and Profs. Joseph Barranco (SF State Univ.), Andisheh Mahdavi (SF State Univ.), and Stephen Kane (UC Riverside).  We also appreciate the assistance of Josh Lamstein, Shervin Sahba, Dirk Kessler, Paul Seawell, Craig Schuler, and Arjun Savel, who helped with various aspects of the code development.  We are grateful to Prof. Eliza Kempton (Univ. of Maryland) for her guidance on installation and usage of the Exo-Transmit code.  We would like to thank Dr. Tom Greene (NASA Ames) for his insight and guidance on a number of issues.  We would also like to acknowledge the assistance of Dr. Natasha Batalha (UC Santa Cruz) the lead author of the PandExo code.
