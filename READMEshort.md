# *JWST* Exoplanet Targeting program (JET)
Computer tool to rank lists of exoplanet targets for atmospheric characterization by *JWST*.

## Description
The *James Webb Space Telescope* (*JWST*) will devote significant observing time to the study of exoplanets. It will not be serviceable as was the Hubble Space Telescope, and therefore the spacecraft/instruments will have a relatively limited life.  It is important to get as much science as possible out of this limited observing time.  We provide a computer tool that can be used to optimize lists of exoplanet targets for atmospheric characterization. Our tools take catalogs of planet detections, either simulated, or actual; categorize the targets by radius and equilibrium temperature; estimate planet masses; generate model spectra and simulated instrument spectra; perform a statistical analysis to determine if the instrument spectra can confirm an atmospheric detection; and finally, rank the targets within each category by observation time required for detection.

## Getting Started

### Hardware Configuration
* The JET code was developed and tested on a modest custom-built desktop computer (see Table 1).  In addition, all production runs to date have been executed on this machine.  

 Table 1. Hardware for Development and Testing
 
	| Component         | Details                |
	|-------------------|------------------------|
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
