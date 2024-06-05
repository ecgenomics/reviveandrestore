#!/bin/bash -ue
/Users/fabio/newcloud/reviveandrestore/mapping.capture/software/bedtools2/bin/intersectBed -abam  test_R_rmdups.bam.stats -b input/tuf.bed > test_R_ontarget.bam
