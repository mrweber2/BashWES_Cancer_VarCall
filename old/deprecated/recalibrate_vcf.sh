#!/bin/bash
#
# recalibrate_vcf.sh <runfile> <sample> <log.in> <log.ou> <qsub>
# 
##redmine=hpcbio-redmine@igb.illinois.edu
redmine=grendon@illinois.edu
if [ $# != 5 ]
then
        MSG="Parameter mismatch. Rerun as: $0 <runfile> <inputfile> <log.in> <log.ou> <qsub> "
        echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        exit 1;
fi
set +x
echo -e "\n\n#####################################################################################"        
echo -e "#############             BEGIN ANALYSIS PROCEDURE                    ###############"
echo -e "#####################################################################################\n\n"        

echo -e "\n\n#####################################################################################"        
echo -e "#############             DECLARING VARIABLES                         ###############"
echo -e "#####################################################################################\n\n"        

set -x
echo `date`
scriptfile=$0
runfile=$1
input_vcf=$2
elog=$3
olog=$4
qsubfile=$5
LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"


if [ ! -s $runfile ]
then
    MSG="$runfile configuration file not found"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "Variant Calling Workflow failure message" "$redmine"
    exit 1;
fi

reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
tmpdir=$( cat $runfile | grep -w TMPDIR | cut -d '=' -f2 )
thr=$( cat $runfile | grep -w PBSCORES | cut -d '=' -f2 )
refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
refgenome=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
dbSNP=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
hapmap=$( cat $runfile | grep -w HAPMAP | cut -d '=' -f2 )
omni=$( cat $runfile | grep -w OMNI | cut -d '=' -f2 )
indels=$( cat $runfile | grep -w INDELS | cut -d '=' -f2 )
phase1=$( cat $runfile | grep -w PHASE1 | cut -d '=' -f2 )
samtools_mod=$( cat $runfile | grep -w SAMTOOLSMODULE | cut -d '=' -f2 )
vcftools_mod=$( cat $runfile | grep -w VCFTOOLSMODULE | cut -d '=' -f2 )
gatkdir=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
ref_local=${refdir}/$refgenome
dbsnp_local=${refdir}/$dbSNP
hapmap_local=${refdir}/$hapmap
omni_local=${refdir}/$omni
indels_local=${refdir}/$indels
phase1_local=${refdir}/$phase1
outputdir=$rootdir
sample=$( basename $input_vcf )
prefix=${sample.GATKCombineGVCF.raw.vcf}
recal_file=${prefix}_recal
tranches_file=${prefix}_tranches
Rplots_snp=${prefix}_plots_snp.R
Rplots_indel=${prefix}_plots_indel.R
outputfile_snp_gvcf=${prefix}_recalibrated_snp_GVCF.vcf
outputfile_indel_gvcf=${prefix}_recalibrated_indel_GVCF.vcf
outputfile_snp_vcf=${prefix}_recalibrated_snp.vcf
outputfile_indel_vcf=${prefix}_recalibrated_indel.vcf
set +x
echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"          	
echo -e "#######   we will need these guys throughout, let's take care of them now   ######"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n"          


echo -e "\n\n##################################################################################" 
echo -e "##################################################################################"        
echo -e "#############                       SANITY CHECK                   ###############"
echo -e "##################################################################################"
echo -e "##################################################################################\n\n"
set -x 
if [ ! -d $tmpdir ]
then
    mkdir -p $tmpdir
fi
if [ ! -d $rootdir ]
then
    MSG="Invalid value specified for OUTPUTDIR=$rootdir in the configuration file."
    echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

if [ `expr ${#input_vcf}` -lt 1 ]
then
    MSG="$input_vcf sample undefined variable"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
    exit 1
fi

if [ ! -s $input_vcf ]
then
    MSG="$input_vcf input file not found"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
    exit 1
fi
if [ `expr ${#sample}` -lt 1 ]
then
    MSG="$input_vcf parsing failed"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
    exit 1
fi
if [ `expr ${#prefix}` -lt 1 ]
then
    MSG="$input_vcf parsing failed"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
    exit 1
fi

set +x
echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"          	
echo -e "#######   STEP1: run VariantRecalibrator on $input_vcf                      ######"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n" 
set -x 
cd $outputdir
set +x
echo -e "\n\n##################################################################################" 
echo -e "########### command one: VariantRecalibrator for SNPs                       #####"
echo -e "##################################################################################\n\n"
set -x

$javadir/java -Xmx8g  -Djava.io.tmpdir=$tmpdir -jar $gatkdir/GenomeAnalysisTK.jar \
         -T VariantRecalibrator \
	 -R $ref_local \
         -input $input_vcf \
         -nt $thr \
	 -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap_local \
         -resource:omni,known=false,training=true,truth=true,prior=12.0   $omni_local \
         -resource:1000G,known=false,training=true,truth=false,prior=10.0 $phase1_local \
	 -resource:dbsnp,known=true,training=false,truth=false,prior=2.0  $dbsnp_local \
         -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff \
	 -mode SNP \
         -recalFile $recal_file \
         -tranchesFile $tranches_file \
         -rscriptFile $Rplots_snp
	 
exitcode=$?
set +x
echo -e "\n\n##################################################################################" 
echo -e "########### command two: sanity check  VariantRecalibrator for SNPs              #####"
echo -e "##################################################################################\n\n"
set -x 
echo `date`
if [ $exitcode -ne 0 ]
then
	 MSG="VariantRecalibrator SNP command failed exitcode=$exitcode. merge for sample $SampleName stopped"
	 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	 exit $exitcode;
fi

if [ ! -s $recal_file ]
then     
    MSG="VariantRecalibrator SNP command did not produce a file $recal_file"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"        
    exit 1;          
fi	

if [ ! -s $tranches_file ]
then     
    MSG="VariantRecalibrator SNP command did not produce a file $tranches_file"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"        
    exit 1;          
fi
set +x
echo -e "\n\n##################################################################################" 
echo -e "########### command three: VariantRecalibrator for Indels                       #####"
echo -e "##################################################################################\n\n"
set -x 

$javadir/java -Xmx8g  -Djava.io.tmpdir=$tmpdir -jar $gatkdir/GenomeAnalysisTK.jar \
         -T VariantRecalibrator \
	 -R $ref_local \
         -input $input_vcf \
         -nt $thr \
         --maxGaussians 4 \
         -resource:mills,known=false,training=true,truth=true,prior=12.0 $mills \
	 -resource:dbsnp,known=true,training=false,truth=false,prior=2.0  $dbsnp_local \
         -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum  -an InbreedingCoeff \
         -mode INDEL \
         -recalFile $recal_file \
         -tranchesFile $tranches_file \
         -rscriptFile $Rplots_indel	 
exitcode=$?
set +x
echo -e "\n\n##################################################################################" 
echo -e "########### command four: sanity check VariantRecalibrator for Indels            #####"
echo -e "##################################################################################\n\n"
set -x 
echo `date`
if [ $exitcode -ne 0 ]
then
	 MSG="VariantRecalibrator INDEL command failed exitcode=$exitcode. merge for sample $SampleName stopped"
	 echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	 exit $exitcode;
fi

if [ ! -s $recal_file ]
then     
    MSG="VariantRecalibrator INDEL command did not produce a file $recal_file"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"        
    exit 1;          
fi	

if [ ! -s $tranches_file ]
then     
    MSG="VariantRecalibrator INDEL command did not produce a file $tranches_file"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"        
    exit 1;          
fi

set +x
echo -e "\n\n##################################################################################"  
echo -e "##################################################################################"          	
echo -e "#######   STEP2: run ApplyRecalibration on $input_vcf                       ######"
echo -e "##################################################################################"  
echo -e "##################################################################################\n\n" 

echo -e "\n\n##################################################################################" 
echo -e "########### command one: ApplyRecalibration for SNPs                       #####"
echo -e "##################################################################################\n\n"
set -x 

$javadir/java -Xmx8g  -Djava.io.tmpdir=$tmpdir -jar $gatkdir/GenomeAnalysisTK.jar \
         -T ApplyRecalibration \
	 -R $ref_local \
         -input $input_vcf \
         -nt $thr \
         --ts_filter_level 99.5 \
         -mode SNP \
         -recalFile $recal_file \
         -tranchesFile $tranches_file \
         -o $outputfile_snp_gvcf
	 
exitcode=$?

set +x
echo -e "\n\n##################################################################################" 
echo -e "########### command three: ApplyRecalibration for INDELs                       #####"
echo -e "##################################################################################\n\n"
set -x

$javadir/java -Xmx8g  -Djava.io.tmpdir=$tmpdir -jar $gatkdir/GenomeAnalysisTK.jar \
         -T ApplyRecalibration \
	 -R $ref_local \
         -input $input_vcf \
         -nt $thr \
         --ts_filter_level 99.0 \
         -mode INDEL \
         -recalFile $recal_file \
         -tranchesFile $tranches_file \
         -o $outputfile_indel_gvcf
	 
exitcode=$?


set +x
echo `date`
echo -e "\n\n##################################################################################"
echo -e "#############    DONE PROCESSING SAMPLE $SampleName. EXITING NOW.  ###############"
echo -e "##################################################################################\n\n"
set -x 

