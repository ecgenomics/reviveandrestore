#!/usr/bin/env nextflow

/*
  _____            _                      _____           _                 
 |  __ \          (_)             ___    |  __ \         | |                
 | |__) |_____   _____   _____   ( _ )   | |__) |___  ___| |_ ___  _ __ ___ 
 |  _  // _ \ \ / / \ \ / / _ \  / _ \/\ |  _  // _ \/ __| __/ _ \| '__/ _ \
 | | \ \  __/\ V /| |\ V /  __/ | (_>  < | | \ \  __/\__ \ || (_) | | |  __/
 |_|  \_\___| \_/ |_| \_/ \___|  \___/\/ |_|  \_\___||___/\__\___/|_|  \___|
                                                                            
                                                                            
  Module name: main.nf
                                                                           
  Endogenous or hDNA Quantification to prepare equi-endogenous pools for capture 
  Based on https://github.com/claudefa/Mapping_Capture_Experiments
*/


// Syntax
nextflow.enable.dsl=2

// Params
params.str = 'Hello world!'
params.logDir = "${baseDir}/log"

println baseDir
println params.logDir

// Processes
process splitLetters {
    output:
    path 'chunk_*'

    """
    printf '${params.str}' | split -b 6 - chunk_
    """
}

process convertToUpper {

    input:
    path x

    output:
    stdout

    """
    cat $x | tr '[a-z]' '[A-Z]'
    """
}

workflow {
    splitLetters | flatten | convertToUpper | view { it.trim() }
}