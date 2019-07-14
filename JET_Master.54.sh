#!/bin/bash

# @author: Charles Fortenbach

# Code associated with paper submitted to AAS Journals (6/10/19)

# Title: A FRAMEWORK FOR OPTIMIZING EXOPLANET TARGET SELECTION FOR
# THE JAMES WEBB SPACE TELESCOPE

# Subprogram:  JET_Master.sh

# Version:  1.0

# This project is licensed under the GNU GPLv3 License.


########## DEFINE CURRENT WORKING DIRECTORY ################
USRDIR=$(pwd)
echo 'USRDIR=$(pwd)' >>~/.bash_profile

########## BLOCK TO SET STELLAR DATA ################
echo 'export PYSYN_CDBS="$USRDIR/grp/hst/cdbs"' >>~/.bash_profile

############ BLOCK TO SET PANDEIA REFERENCE DATA #########################
echo 'export pandeia_refdata="$USRDIR/pandeia_data-1.4"' >>~/.bash_profile

############## MAKE SURE YOUR BASH PROFILE IS SOURCED ################
source ~/.bash_profile

########## RUN ###############
python ExoT_Master.57.py
python Pdxo_Master.66.py
python Rank_Master.08.py
