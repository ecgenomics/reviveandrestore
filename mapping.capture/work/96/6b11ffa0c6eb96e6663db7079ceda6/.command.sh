#!/bin/bash -ue
/Users/fabio/newcloud/new/mapping.capture/software/samtools-1.20/samtools sort -o test_R_rmdups.sorted.bam test_R_rmdups.bam; 	/Users/fabio/newcloud/new/mapping.capture/software/samtools-1.20/samtools index test_R_rmdups.bam;
/Users/fabio/newcloud/new/mapping.capture/software/samtools-1.20/samtools view -h -q 30 -F 256 test_R_rmdups.sorted.bam | samtools view -hb > test_R_rmdups.qual.bam; 	/Users/fabio/newcloud/new/mapping.capture/software/samtools-1.20/samtools index test_R_rmdups.qual.bam;
