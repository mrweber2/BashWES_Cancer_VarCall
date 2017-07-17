#!/bin/bash
#
# index_bams.sh
# samtools index recalibrated bams, we need them to run MuTect2

set -x
echo `date`

allbams=$1

if [ ! -s $allbams ]
then
     echo "$allbams file not found"
     exit 1
fi

module load samtools/1.3

while read bam
do
   echo -e "Processign next bam=$bam"
   
   dirname=$( dirname $bam )
   cd $dirname
   samtools index $bam

done < $allbams