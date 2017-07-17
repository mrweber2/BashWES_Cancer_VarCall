#!/bin/bash

scriptfile=$0
runfile=$1
SampleName=$2
chr=$3
elog=$4
olog=$5
qsubfile=$6
LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
email=$( cat $runfile | grep -w EMAIL | cut -d '=' -f2 )
rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
tmpdir=$( cat $runfile | grep -w TMPDIR | cut -d '=' -f2 ) 
sampleinfo=$( cat $runfile | grep SAMPLEINFORMATION | cut -d '=' -f2)
refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
refgenome=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
dbSNP=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 ) 
javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
gatkdir=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
samtoolsdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
rdir=$( cat $runfile | grep -w RDIR | cut -d '=' -f2 )

thr=$( cat $runfile | grep -w PBSCORES | cut -d '=' -f2 )

ref_local=${refdir}/$refgenome
dbsnp_local=${refdir}/$dbSNP
indel_local=${refdir}/$indeldir
echo -e "\n\n##################################################################################" >&2 
echo -e "##################################################################################" >&2    
echo -e "#######   we will need these guys throughout, let's take care of them now   ######" >&2
echo -e "##################################################################################" >&2 
echo -e "##################################################################################\n\n" >&2

outputdir=$rootdir/$SampleName
AlignDir=$outputdir/align
RealignDir=$outputdir/realign
inputbam=$AlignDir/${SampleName}.wdups.sorted.bam                # name of the bam file that align-dedup produced
dedupsortedbam=${SampleName}.${chr}.wdups.sorted.bam               # name of the aligned file -one chr out of inputbam



echo -e "\n\n##################################################################################"


echo -e "\n\n########################################################################################"
echo -e "##############                 plotting now for this chr=$chr                            ##################"	
echo -e "########################################################################################\n\n"

module load R/3.2.0
Rscript ${scriptdir}/ParameterSweep/bqsrplotting.R ${rootdir}/sweepBQSR.${SampleName}/ $thr $chr



echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
