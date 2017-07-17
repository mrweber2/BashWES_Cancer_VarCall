#!/bin/bash
#PBS -q default
#PBS -m ae
#PBS -M grendon@illinois.edu
#PBS -l nodes=1:ppn=24
#PBS -N MuTec_Adjacent_Normal_Sample_A_vs_Cancer_1_Sample_A
#PBS -o /home/groups/hpcbio/projects/LifeSpan/exome-March2016/src/qsub.MuTec_Adjacent_Normal_Sample_A_vs_Cancer_1_Sample_A.ou
#PBS -e /home/groups/hpcbio/projects/LifeSpan/exome-March2016/src/qsub.MuTec_Adjacent_Normal_Sample_A_vs_Cancer_1_Sample_A.er


set -x
echo `date`

set +x
echo "\n\n##############  Step1: Setup\n\n"
set -x

refgenome=/home/groups/hpcbio/projects/LifeSpan/exome-March2016/genome/genome.fa
dbSNP=/home/groups/hpcbio/projects/LifeSpan/exome-March2016/genome/dbSNP.vcf
normal_vcf=/home/groups/hpcbio/projects/LifeSpan/exome-March2016/results/2016-05-12-exome2-wpicard/delivery/Exome2_Adjacent_Normal_Sample_A/Exome2_Adjacent_Normal_Sample_A.sorted.raw.vcf
normal=/home/groups/hpcbio/projects/LifeSpan/exome-March2016/results/2016-05-12-exome2-wpicard/delivery/Exome2_Adjacent_Normal_Sample_A/Exome2_Adjacent_Normal_Sample_A.recalibrated.bam
tumor=/home/groups/hpcbio/projects/LifeSpan/exome-March2016/results/2016-05-12-exome2-wpicard/delivery/Exome2_Cancer_1_Sample_A/Exome2_Cancer_1_Sample_A.recalibrated.bam
outputdir=/home/groups/hpcbio/projects/LifeSpan/exome-March2016/results/2016-05-16-MuTec-pilotrun
outputfile=MuTect2.Adjacent_Normal_Sample_A_vs_Cancer_1_Sample_A.vcf
gatk_dir=/home/apps/gatk/gatk-3.5
tmpdir=/state/partition1/grendon
thr=23
module load gatk/3.5

if [ ! -d $outputdir ]
then
    mkdir -p $outputdir
fi

cd $outputdir 

exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
    MSG="Invalid value specified for outputdir=$outputdir  exiting now."
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" #| mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;

fi
set +x
echo "\n\n##############  Step2: MuTect command with -PON param \n\n"
set -x

#java -Xmx10g  -Djava.io.tmpdir=$tmpdir -jar $gatk_dir/GenomeAnalysisTK.jar \
#         -T MuTect2 \
#	 -R $refgenome \
#	 --dbsnp $dbSNP \
#	 -I:tumor $tumor \
#	 -I:normal $normal \
#	 --normal_panel $normal_vcf \
#	 -nct $thr \
#	 -o $outputfile

set +x
echo "\n\n##############  Step2: MuTect command without -PON param \n\n"
set -x

java -Xmx10g  -Djava.io.tmpdir=$tmpdir -jar $gatk_dir/GenomeAnalysisTK.jar \
         -T MuTect2 \
	 -R $refgenome \
	 --dbsnp $dbSNP \
	 -I:tumor $tumor \
	 -I:normal $normal \
	 -nct $thr \
	 -o $outputfile
	 
exitcode=$?
set +x
echo "\n\n##############  Step3: Sanity Check \n\n"
set -x 

echo `date`
if [ $exitcode -ne 0 ]
then
	 MSG="MuTect2 command failed exitcode=$exitcode\nInputs:\n\nnormal=$normal\tumor=$tumor\n"
	 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
	 exit $exitcode;
fi

if [ ! -s $outputfile ]
then
	 MSG="MuTect2 produced an empty file exitcode=$exitcode\nInputs:\n\nnormal=$normal\tumor=$tumor\n"
	 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
	 exit $exitcode;
else 
set +x
echo "\n\n##############  Done. Exiting now \n\n"
set -x
