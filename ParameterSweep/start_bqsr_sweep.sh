#!/bin/bash

################################################################################################ 
# Program to do a parameter sweep from human samples of WES short reads
# In order to run this pipeline please type at the command line
# start_bqsr_sweep.sh <runfile>
################################################################################################

set -x
redmine=hpcbio-redmine@igb.illinois.edu
scriptfile=$0
runfile=$1

if [ $# != 1 ] 
then
     MSG="Parameter mismatch.\nRerun like this: $0 <runfile>\n"
     echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "BQSR parameter sweep Workflow failure message" "$redmine"
     exit 1;
fi

reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
tmpdir=$( cat $runfile | grep -w TMPDIR | cut -d '=' -f2 )
deliverydir=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 )
scriptdir=$( cat $runfile | grep -w SCRIPTDIR | cut -d '=' -f2 )
email=$( cat $runfile | grep -w EMAIL | cut -d '=' -f2 )
sampleinfo=$( cat $runfile | grep -w SAMPLEINFORMATION | cut -d '=' -f2 )
numsamples=$(wc -l $sampleinfo)
refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
refgenome=$( cat $runfile | grep -w REFGENOME | cut -d '=' -f2 )
dbSNP=$( cat $runfile | grep -w DBSNP | cut -d '=' -f2 )
sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )
dup_cutoff=$( cat $runfile | grep -w  DUP_CUTOFF | cut -d '=' -f2 )
map_cutoff=$( cat $runfile | grep -w  MAP_CUTOFF | cut -d '=' -f2 )
indices=$( cat $runfile | grep -w CHRNAMES | cut -d '=' -f2 | tr ':' ' ' )
analysis=$( cat $runfile | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
alignertool=$( cat $runfile | grep -w ALIGNERTOOL | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
markduplicates=$( cat $runfile | grep -w MARKDUPLICATESTOOL | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
samblasterdir=$( cat $runfile | grep -w SAMBLASTERDIR | cut -d '=' -f2 )
picardir=$( cat $runfile | grep -w PICARDIR | cut -d '=' -f2 )
gatkdir=$( cat $runfile | grep -w GATKDIR | cut -d '=' -f2 )
samtoolsdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
bwamemdir=$( cat $runfile | grep -w BWAMEMDIR | cut -d '=' -f2 )
javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
novocraftdir=$( cat $runfile | grep -w NOVOCRAFTDIR | cut -d '=' -f2 )
fastqcdir=$( cat $runfile | grep -w FASTQCDIR | cut -d '=' -f2 )
thr=$( cat $runfile | grep -w PBSCORES | cut -d '=' -f2 )
nodes=$( cat $runfile | grep -w PBSNODES | cut -d '=' -f2 )
queue=$( cat $runfile | grep -w PBSQUEUE | cut -d '=' -f2 )
pbswalltime=$( cat $runfile | grep -w PBSWALLTIME | cut -d '=' -f2 )

TopOutputLogs=$outputdir/logs

generic_qsub_header=$TopOutputLogs/qsubGenericHeader
truncate -s 0 $generic_qsub_header
echo "#!/bin/bash" > $generic_qsub_header
echo "#PBS -q $queue" >> $generic_qsub_header
echo "#PBS -m ae" >> $generic_qsub_header
echo "#PBS -M $email" >> $generic_qsub_header
echo "#PBS -l nodes=$nodes:ppn=$thr" >> $generic_qsub_header
echo "#PBS -l walltime=${pbswalltime}" >> $generic_qsub_header

while read sampleLine 
do
    if [ `expr ${#sampleLine}` -lt 1 ] ;  then
        set +x 
        echo -e "\n\n########################################################################################" >&2
        echo -e "##############                 skipping empty line        ##############################" >&2
        echo -e "########################################################################################\n\n" >&2
    fi
    sample=$( echo "$sampleLine" | cut -d ' ' -f 1 )
    for chr in $indices
    do
        set +x
        echo -e "\n\n########################################################################################" >&2 
        echo -e "####   BQSR sweep script for SAMPLE $sample chr=$chr                              #######" >&2
        echo -e "########################################################################################\n\n" >&2
        set -x
        qsub1=$TopOutputLogs/qsub.bqsr_sweep.$sample.$chr
        cat $generic_qsub_header > $qsub1
        echo "#PBS -N bqsr_sweep.$sample.$chr" >> $qsub1
        echo "#PBS -o $TopOutputLogs/log.bqsr_sweep.$sample.$chr.ou" >> $qsub1
        echo "#PBS -e $TopOutputLogs/log.bqsr_sweep.$sample.$chr.in" >> $qsub1
        
        echo "$scriptdir/ParameterSweep/bqsr_sweep.sh $runfile $sample $chr $TopOutputLogs/log.bqsr_sweep.$sample.$chr.in $TopOutputLogs/log.bqsr_sweep.$sample.$chr.ou $TopOutputLogs/qsub.bqsr_sweep.$sample.$chr" >> $qsub1
        `chmod a+r $qsub1`               
        bqsr_sweepjobid=`qsub $qsub1` 
        echo $bqsr_sweepjobid >> $TopOutputLogs/pbs.BQSR_SWEEP.$sample
        echo $bqsr_sweepjobid >> $TopOutputLogs/pbs.summary_dependencies
        
        echo `date`
        if [ `expr ${#bqsr_sweepjobid}` -lt 1 ]
        then
                MSG="unable to launch qsub bqsr_sweep job for $sample. Exiting now"
                echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
                exit 1        
        fi
    done
    bqsrjobids=$( cat $TopOutputLogs/pbs.BQSR_SWEEP.$sample | sed "s/\.[a-z]*//g" | tr "\n" ":" )

done <  $sampleinfo

echo -e "########################################################################################" >&2
echo -e "#################     Now, we need to generate summary           #######################" >&2
echo -e "########################################################################################" >&2
echo -e "########################################################################################\n\n" >&2
    
: << 'end_comment_out'
alljobids=$( cat $TopOutputLogs/pbs.summary_dependencies | sed "s/\.[a-z]*//g" | tr "\n" ":" )

set +x
echo -e "\n\n### this list of jobids=[$alljobids] will be used to hold execution of summary.sh #####\n\n" >&2
set -x

qsub2=$TopOutputLogs/qsub.summary
cat $generic_qsub_header > $qsub2
echo "#PBS -N Summary_vcall" >> $qsub2
echo "#PBS -o $TopOutputLogs/log.summary.ou" >> $qsub2
echo "#PBS -e $TopOutputLogs/log.summary.in" >> $qsub2
echo "#PBS -W depend=afterok:$alljobids " >> $qsub2
echo "$scriptdir/summary.sh $runfile $TopOutputLogs/log.summary.in $TopOutputLogs/log.summary.ou $TopOutputLogs/qsub.summary" >> $qsub2

`chmod a+r $qsub2`
lastjobid=`qsub $qsub2`
echo $lastjobid >> $TopOutputLogs/pbs.SUMMARY
echo `date`


if [ `expr ${#lastjobid}` -lt 1 ]
then
     MSG="unable to launch qsub summary job. Exiting now"
     echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
     exit 1
fi

end_comment_out
