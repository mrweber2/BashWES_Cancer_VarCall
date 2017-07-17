#!/bin/bash
redmine=hpcbio-redmine@igb.illinois.edu
##redmine=grendon@illinois.edu
if [ $# != 4 ]
then
        MSG="Parameter mismatch."
        echo -e "program=$0 stopped. Reason=$MSG" | mail -s 'Variant Calling Workflow failure message' "$redmine"
        
        exit 1;
fi


set -x
echo `date`	
scriptfile=$0
runfile=$1
elog=$2
olog=$3
qsubfile=$4
LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
sampleinfo=$( cat $runfile | grep -w SAMPLEINFORMATION | cut -d '=' -f2 )        
reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
deliverydir=$( cat $runfile | grep -w DELIVERYFOLDER | cut -d '=' -f2 )
summaryok="YES"
samplesProcessed=""
set +x
echo -e "\n\n############################################################################################################"
echo -e "###########                        SANITY   CHECK   HERE                                  ##################"
echo -e "############################################################################################################\n\n"
set -x 
if [ ! -s $outputdir/$deliverydir/docs/Summary.Report ]
then 
    MSG="$outputdir/$deliverydir/docs/Summary.Report  file not found"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"                     
    exit 1                
fi

while read sampleLine
do
    if [ `expr ${#sampleLine}` -lt 1 ]
    then
	set +x
	echo -e "\n\n############################################################################################################"
	echo -e "##############                              skipping empty line        #####################################"
	echo -e "############################################################################################################\n\n"
	set -x 
    else
	set +x
	echo -e "\n\n############################################################################################################"
	echo -e "##############                      Processing next line $sampleLine   #####################################"
	echo -e "############################################################################################################\n\n"
	set -x 
	sample=$( echo "$sampleLine" | cut -d ' ' -f 1 ) 
	set +x
	echo -e "\n\n############################################################################################################"
	echo -e "############## checking that there are results available for sample $sample in $outputdir/$sample/ #########"
	echo -e "############################################################################################################\n\n"
	set -x 
	if [ -s $outputdir/$sample ]
	then
	    echo -e "$sample\tResults folder:\t$outputdir/$sample" >> $samplesProcessed
	    echo -e "########## All results seem to be OK $outputdir/$sample"
	else
	    echo -e "$sample\tResults folder:\tEXECUTION FAILED TO PRODUCE SOME RESULTS" >> $samplesProcessed 
	    echo -e "########## Some results seem to be MISSING $outputdir/$sample"                    
	    summaryok="NO"
	fi
   fi 
done <  $sampleinfo                

set +x
echo -e "\n\n############################################################################################################"
echo -e "####   now putting together the second part of the Summary.Report file with the list of jobs executed   ####"
echo -e "############################################################################################################\n\n"
set -x 
listjobids=$( cat $outputdir/logs/pbs.* | sort | uniq | tr "\n" "\t" )

if [ $summaryok == "YES" ]
then
    MSG="Variant calling workflow run by username: $USER finished with ALL  jobs with exit code 0 at: "$( echo `date` )	
else 
    MSG="Variant calling workflow run by username: $USER finished with SOME  jobs with exit code 0 at: "$( echo `date` )
fi

LOGS="Results and execution logs can be found at $outputdir\n\nDeliverables are located at $outputdir/$deliverydir\n\nJOBIDS\n\n$listjobids\n\nThis jobid:${PBS_JOBID}\n\n"
echo -e "$MSG\n\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
echo -e "$MSG\n\n$LOGS" >> $outputdir/$deliverydir/docs/Summary.Report

chmod -R g+r $outputdir
set +x
echo -e "\n\n############################################################################################################"
echo -e "############                              DONE. EXITING NOW                                      ###########"
echo -e "############################################################################################################\n\n"
set -x
