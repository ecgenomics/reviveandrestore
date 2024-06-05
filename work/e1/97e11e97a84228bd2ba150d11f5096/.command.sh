#!/bin/bash -ue
/Users/fabio/newcloud/new/software/bwa/bwa mem -t 4 -M input/hg19.fa.masked pL4811A1_hP37_trimmed_R1.fastq pL4811A1_hP37_trimmed_R2.fastq | 
awk 'substr($0,1,1)=="@" || ($9>= 35 && $9<=500) || ($9<=-35 && $9>=-500)' |
/Users/fabio/newcloud/new/software/samtools-1.20/samtools view -Shb - |
java -Djava.io.tmpdir=/Users/fabio/newcloud/new/output/temp -jar /Users/fabio/newcloud/new/software/picard/build/libs/picard.jar AddOrReplaceReadGroups I=/dev/stdin O=pL4811A1_hP37_pe.bam     SORT_ORDER=coordinate QUIET=TRUE COMPRESSION_LEVEL=9 MAX_RECORDS_IN_RAM=150000     RGLB=pL4811A1_hP37 RGPL=Illumina RGPU=pL4811A1_hP37 RGSM=pL4811A1_hP37 RGID=pL4811A1_hP37    VALIDATION_STRINGENCY=SILENT
