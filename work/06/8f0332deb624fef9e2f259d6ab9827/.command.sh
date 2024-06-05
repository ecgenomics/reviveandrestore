#!/bin/bash -ue
/Users/fabio/newcloud/new/mapping.capture/software/bwa/bwa mem -t 4 -M input/hg19.fa.masked test_R_trimmed_R1.fastq test_R_trimmed_R2.fastq | 
awk 'substr($0,1,1)=="@" || ($9>= 35 && $9<=500) || ($9<=-35 && $9>=-500)' |
/Users/fabio/newcloud/new/mapping.capture/software/samtools-1.20/samtools view -Shb - |
java -Djava.io.tmpdir=/Users/fabio/newcloud/new/mapping.capture/output/temp -jar /Users/fabio/newcloud/new/mapping.capture/software/picard/build/libs/picard.jar AddOrReplaceReadGroups I=/dev/stdin O=test_R_pe.bam     SORT_ORDER=coordinate QUIET=TRUE COMPRESSION_LEVEL=9 MAX_RECORDS_IN_RAM=150000     RGLB=test_R RGPL=Illumina RGPU=test_R RGSM=test_R RGID=test_R    VALIDATION_STRINGENCY=SILENT
