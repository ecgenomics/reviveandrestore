#!/usr/bin/env python3

#  _____            _                      _____           _                 
# |  __ \          (_)             ___    |  __ \         | |                
# | |__) |_____   _____   _____   ( _ )   | |__) |___  ___| |_ ___  _ __ ___ 
# |  _  // _ \ \ / / \ \ / / _ \  / _ \/\ |  _  // _ \/ __| __/ _ \| '__/ _ \
# | | \ \  __/\ V /| |\ V /  __/ | (_>  < | | \ \  __/\__ \ || (_) | | |  __/
# |_|  \_\___| \_/ |_| \_/ \___|  \___/\/ |_|  \_\___||___/\__\___/|_|  \___|
#                                                                            
#                                                                            
#  Script name: cfolders.py
#                                                                           
#  This creates the folder structure for all the analyses. It is called by
#  main.nf and imports the global settings.

import sys, shutil
from pathlib import Path as p

session = sys.argv[1]

folderlist = []

# DEFINE FOLDERS: Demultiplexed fastq folders
#######################################################

out = session + "/fastqs/demultiplex/out"                                   # Output fastq demultiplex
qu = session + "/fastqs/demultiplex/qu"                                     # output fastq qu

# Add to folderlist
folderlist.append(out)
folderlist.append(qu)

# DEFINE FOLDERS: Trimmed adapters outputs
#######################################################

tr_out = session + "/fastqs/trimmed/demultiplex/out"                        # Output trimmed out 
tr_qu = session + "/fastqs/trimmed/demultiplex/qu"                          # Output trimmed qu

# Add to folderlist
folderlist.append(tr_out)
folderlist.append(tr_qu)

# DEFINE FOLDERS: Mapping HG19 hDNA
#######################################################

bam_hg19 = session + "/bam_hg19"                                              # Bam HG19 output

# Add to folderlist
folderlist.append(bam_hg19)

# DEFINE FOLDERS: Remove duplicates
#######################################################

bam_hg19_no_dups = session + "/bam_hg19_no_dups"

# Add to folderlist
folderlist.append(bam_hg19_no_dups)

# OVERWRITE SESSION FOLDER
#######################################################

if p(session).is_dir():                                                     ### DEV!!!
    shutil.rmtree(session)

# CREATE FOLDERS: all folders
######################################################

for x in folderlist:
    p(x).mkdir(parents=True)