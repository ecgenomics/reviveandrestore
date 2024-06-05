#!/bin/bash -ue
java -Xmx8g -Djava.io.tmpdir=/Users/fabio/newcloud/new/mapping.capture/output/temp -jar /Users/fabio/newcloud/new/mapping.capture/software/picard/build/libs/picard.jar MarkDuplicates     I=test_R_srtd.bam     O=test_R_rmdups.bam     MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 METRICS_FILE=test_R_rmdups.bam.stats      REMOVE_DUPLICATES=true ASSUME_SORTED=true MAX_RECORDS_IN_RAM=1500000     VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true COMPRESSION_LEVEL=9
