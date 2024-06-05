#!/bin/bash -ue
/Users/fabio/newcloud/reviveandrestore/mapping.capture/software/samtools-1.20/samtools sort -T test_R -o test_R_rmdups.sorted.bam; 	/Users/fabio/newcloud/reviveandrestore/mapping.capture/software/samtools-1.20/samtools index test_R_rmdups.bam;
/Users/fabio/newcloud/reviveandrestore/mapping.capture/software/samtools-1.20/samtools view -h -q 30 -F 256 -F test_R_rmdups.sorted.bam | samtools view -hb > test_R_rmdups.qual.bam; 	/Users/fabio/newcloud/reviveandrestore/mapping.capture/software/samtools-1.20/samtools index test_R_rmdups.qual.bam;
