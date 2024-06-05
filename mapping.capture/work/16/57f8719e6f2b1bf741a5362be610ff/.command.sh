#!/bin/bash -ue
/Users/fabio/newcloud/reviveandrestore/mapping.capture/software/bwa/bwa mem -t 4 -M input/hg19.fa.masked test_R.merged.fastq |
/Users/fabio/newcloud/reviveandrestore/mapping.capture/software/samtools-1.20/samtools view -Sbh - |
java -Djava.io.tmpdir=/Users/fabio/newcloud/reviveandrestore/mapping.capture/output/temp     -jar /Users/fabio/newcloud/reviveandrestore/mapping.capture/software/picard/build/libs/picard.jar AddOrReplaceReadGroups     I=/dev/stdin O=test_R_se.bam     SORT_ORDER=coordinate QUIET=TRUE COMPRESSION_LEVEL=9     MAX_RECORDS_IN_RAM=150000 RGLB=test_R     RGPL=Illumina RGPU=test_R RGSM=test_R     RGID=test_R VALIDATION_STRINGENCY=SILENT
