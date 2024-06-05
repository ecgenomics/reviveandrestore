#!/bin/bash -ue
/Users/fabio/newcloud/new/software/samtools-1.20/samtools sort -o pL4811A1_hP37_rmdups.sorted.bam pL4811A1_hP37_rmdups.bam; 	/Users/fabio/newcloud/new/software/samtools-1.20/samtools index pL4811A1_hP37_rmdups.bam;
/Users/fabio/newcloud/new/software/samtools-1.20/samtools view -h -q 30 -F 256 pL4811A1_hP37_rmdups.sorted.bam | samtools view -hb > pL4811A1_hP37_rmdups.qual.bam; 	/Users/fabio/newcloud/new/software/samtools-1.20/samtools index pL4811A1_hP37_rmdups.qual.bam;
