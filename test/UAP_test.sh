#!/bin/bash

workdir=`pwd`
HG19FA=$1
UAPDIR=$(dirname `pwd`)
export UAP_HOME=$UAPDIR

# AnnoFastqWithUMI
${UAP_HOME}/bin/UAP AnnoFastqWithUMI -f Test_1.fq.gz -r Test_2.fq.gz -o Test -d $workdir/Split_test -e 0

# BamConsensusAnalysis
${UAP_HOME}/bin/UAP BamConsensusAnalysis -i Test.bam -o Test -d $workdir/Bam_test -b test.bed -F $HG19FA -m 10 -s 3 -u Duplex -B 0.1 -r 1

# ErrorRateStats
mkdir -p Error_stat_test
${UAP_HOME}/bin/UAP ErrorRateStats -i Test.bam -o Error_stat_test/Test -F $HG19FA -m 10
