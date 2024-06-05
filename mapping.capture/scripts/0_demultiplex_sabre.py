#!/usr/bin/env python3

#  _____            _                      _____           _                 
# |  __ \          (_)             ___    |  __ \         | |                
# | |__) |_____   _____   _____   ( _ )   | |__) |___  ___| |_ ___  _ __ ___ 
# |  _  // _ \ \ / / \ \ / / _ \  / _ \/\ |  _  // _ \/ __| __/ _ \| '__/ _ \
# | | \ \  __/\ V /| |\ V /  __/ | (_>  < | | \ \  __/\__ \ || (_) | | |  __/
# |_|  \_\___| \_/ |_| \_/ \___|  \___/\/ |_|  \_\___||___/\__\___/|_|  \___|
#                                                                            
#                                                                            
#  Script name: 0_demultiplex_sabre.py (this will change)
#  The first step, that consists in de-multiplexing the fastq files.
#  

import sys, os, glob

### Sabre launch wrapper FUNCTION
def sabre(exec, mismatches, file1, file2, barcode_file, unknown_output1, unknown_output2):
    command = f"{exec} pe -m {mismatches} -f {file1} -r{file2} -b {barcode_file} -u {unknown_output1} -w{unknown_output2}"
    os.system(command)

### Parameters
input_folder  = "../input.provo"   # Bound to change
output_folder = sys.argv[1]

### Processing parameters
inputs = glob.glob(input_folder + "/*_1.fastq.gz")

### Run if main
if __name__ == "__main__":
    exec = ""
    mismatches = ""
    file1 = ""
    file2 = ""
    barcode_file = "",
    unknown_output1 = ""
    unknown_output2 = ""

    sabre(mismatches, file1, file2, barcode_file, unknown_output1, unknown_output2)