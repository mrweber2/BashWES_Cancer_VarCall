#!/bin/bash

set -x
echo `date`

module load bcftools/1.3

cd /home/groups/hpcbio/projects/LifeSpan/exome-March2016/results
mkdir -p 2016-05-14-VCFStats
cd 2016-05-14-VCFStats

truncate -s 0 VCF_Stats.txt

echo -e "###### produce vcf stats with bvftools for all vcfs in delivery folders"

while read raw_vcf
do
    if [ ! -s $raw_vcf ]
    then
        echo -e "\n\n### $raw_vcf file not found. skipping it####\n\n" >> VCF_Stats.txt
    else
        echo -e "\n\n### VCF STATS FOR: $raw_vcf ####\n\n" >> VCF_Stats.txt
        bcftools stats $raw_vcf > ${raw_vcf}.stats

        stats=${raw_vcf}.stats

        if [ ! -s $stats ]
        then
           echo "$stats file not found" >> VCF_Stats.txt
        else
           grep "^SN" $stats >> VCF_Stats.txt
        fi
    fi
    echo `date`
    #sleep 5
done < /home/groups/hpcbio/projects/LifeSpan/exome-March2016/src/all.vcf.filenames

echo `date`
echo "###### DONE. exiting now"

