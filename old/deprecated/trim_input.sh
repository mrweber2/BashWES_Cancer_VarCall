#!/bin/sh

# This script checks the quality of reads in a folder, and trims known adapter sequences from reads using trimmomatic. It should be called as:
# trim_input.sh <runfile>

runfile=$1
redmine=hpcbio-redmine@igb.illinois.edu

if [ $# != 1 ]
then
     MSG="Parameter mismatch.\nRerun like this: $0 <runfile>\n"     
     echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "Variant Calling Workflow failure message" "$redmine"
     exit 1;
else
    set -x
    echo `date`
    rootdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
    tmpdir=$( cat $runfile | grep -w TMPDIR | cut -d '=' -f2 )
    email=$( cat $runfile | grep -w EMAIL | cut -d '=' -f2 )
    reportticket=$( cat $runfile | grep -w REPORTTICKET | cut -d '=' -f2 )
    sampleinfo=$( cat $runfile | grep SAMPLEINFORMATION | cut -d '=' -f2)
    adapters=$( cat $runfile | grep ADAPTERS | cut -d '=' -f2 )
    fastqcdir=$( cat $runfile | grep -w FASTQCDIR | cut -d '=' -f2)
    trimmomaticdir=$( cat $runfile | grep -w TRIMMOMATICDIR | cut -d '=' -f2)
    trimmomaticparams=$( cat $runfile | grep -w TRIMMOMATICPARAMS | cut -d '=' -f2)
    javadir=$( cat $runfile | grep -w JAVADIR | cut -d '=' -f2 )
    TopOutputLogs=$rootdir/logs
    thr=$( cat $runfile | grep -w PBSCORES | cut -d '=' -f2 )
    nodes=$( cat $runfile | grep -w PBSNODES | cut -d '=' -f2 )
    queue=$( cat $runfile | grep -w PBSQUEUE | cut -d '=' -f2 )
    pbswalltime=$( cat $runfile | grep -w PBSWALLTIME | cut -d '=' -f2 )

    echo -e "\n\n############ parameter check                  #################\n\n"
    set -x 
    if [ ! -s $sampleinfo ]
    then
	echo -e "$0 stopped at line $LINENO. \nREASON=input file not found $sampleinfo"| mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
    fi
    if [ ! -s $adapters ]
    then
	echo -e "$0 stopped at line $LINENO. \nREASON=input file not found $adapters"| mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
    fi
    if [ ! -d $rootdir ]; then
        mkdir $rootdir
    else
	    if [ `expr ${#rootdir}` -gt 1 ]  
		    rm -rf $rootdir/*
	    else
		    echo -e "$0 stopped at line $LINENO. \nREASON=ROOTDIR not specified correctly"| mail -s "[Task #${reportticket}]" "$redmine,$email"
	    	    exit 1;
	    fi
    fi

    setfacl -Rm   g::rwx $rootdir  #gives the group rwx permission, and to subdirectories
    setfacl -Rm d:g::rwx $rootdir  #passes the permissions to newly created files/folders

    if [ ! -d $TopOutputLogs ]
    then
	echo -e "creating Logs folder $TopOutputLogs";
	`mkdir $TopOutputLogs`
    fi
    if [ ! -d $tmpdir ]
    then
	echo -e "creating tmp folder $tmpdir";
        `mkdir $tmpdir`
    fi

    set +x
    echo -e "\n\n############ parameters ok                  #################\n\n"
    set -x

    generic_qsub_header=$TopOutputLogs/qsubGenericHeader_qc_trimming
    truncate -s 0 $generic_qsub_header
    echo "#!/bin/bash" > $generic_qsub_header
    echo "#PBS -q $queue" >> $generic_qsub_header
    echo "#PBS -m ae" >> $generic_qsub_header
    echo "#PBS -M $email" >> $generic_qsub_header
    echo "#PBS -l nodes=$nodes:ppn=$thr" >> $generic_qsub_header
    echo "#PBS -l walltime=${pbswalltime}" >> $generic_qsub_header
 
    echo -e "\n\n############ trimming loop start here       #################\n\n"

    while read sampleLine
    do
       set +x	
       echo -e "\n############# processing next line in file...#################\n"
       set -x 
       if [ `expr ${#sampleLine}` -lt 1 ]
       then
	   set +x
           echo -e "\n############ skipping empty line            #################\n"
       else
           echo -e "\n############ processing: $sampleLine     #################\n"
           set -x 
           # we are assuming input is paired reads, three fields per row as follows

           samplename=$( echo "$sampleLine" | cut -d ' ' -f1 )
           R1=$( echo "$sampleLine" | cut -d ' ' -f2 )
           R2=$( echo "$sampleLine" | cut -d ' ' -f3 )
          
	   if [ `expr ${#samplename}` -lt 1 ]
 	   then
	       echo -e "$0 stopped at line $LINENO\nREASON=samplename string not found $samplename"| mail -s "[Task #${reportticket}]" "$redmine,$email"
	       exit 1;
           fi

	   if [ ! -s $R1 ]
	   then
	       echo -e "$0 stopped at line $LINENO\nREASON=reads file not found $R1"| mail -s "[Task #${reportticket}]" "$redmine,$email"
	       exit 1;
           fi

	   if [ ! -s $R2 ]
	   then
	       echo -e "$0 stopped at line $LINENO\nREASON=reads file not found $R2"| mail -s "[Task #${reportticket}]" "$redmine,$email"
	       exit 1;
           fi
	   set +x
           echo -e "\n############ preparing output folders      #################\n"
	   set -x 
           outputdir=$rootdir/${samplename}/trimmed
           fqdir2=$rootdir/${samplename}/FastQC-trimmed
           fqdir1=$rootdir/${samplename}/FastQC-raw
           b1=${samplename}_R1           
           b2=${samplename}_R2


           mkdir -p $fqdir1
           mkdir -p $fqdir2
           mkdir -p $outputdir          
           set +x 
           echo -e "\n############ results will go here: outputdir=$outputdir    #################\n"           
	   set -x 
 

	   qsub1=$TopOutputLogs/qsub.trim.$samplename
	   cat $generic_qsub_header > $qsub1
	   echo "#PBS -V" >> $qsub1
	   echo "#PBS -N trim.${samplename}" >> $qsub1
	   echo "#PBS -e $TopOutputLogs/qsub.trim.${samplename}.er" >> $qsub1
	   echo "#PBS -o $TopOutputLogs/qsub.trim.${samplename}.ou" >> $qsub1
           echo "set -x" >> $qsub1
           echo "echo step1 fastqc on raw reads $R1 $R2" >> $qsub1           

           echo "$fastqcdir/fastqc -o $fqdir1 -t $thr $R1" >> $qsub1
           echo "$fastqcdir/fastqc -o $fqdir1 -t $thr $R2" >> $qsub1
           echo "echo `date`" >>  $qsub1           
           echo "echo step2 trim raw reads" >> $qsub1           
	   echo "$javadir/java -jar $trimmomaticdir PE\
		   -threads $thr \
		   -trimlog $outputdir/${samplename}_trim.log \
		   $R1 $R2 \
		   $outputdir/${b1}.paired.fq.gz $outputdir/${b1}.unpaired.fq.gz \
		   $outputdir/${b2}.paired.fq.gz $outputdir/${b2}.unpaired.fq.gz \
		   ILLUMINACLIP:${adapters}${trimmomaticparams} " >> $qsub1
           echo "echo `date`" >>  $qsub1
           echo "echo step3 fastqc on trimmed reads" >> $qsub1           
           echo "$fastqcdir/fastqc -o $fqdir2 -t $thr $outputdir/${b1}.paired.fq.gz" >> $qsub1
           echo "$fastqcdir/fastqc -o $fqdir2 -t $thr $outputdir/${b2}.paired.fq.gz" >> $qsub1
           echo "echo `date`" >>  $qsub1
           echo "echo exiting now" >> $qsub1           
	   echo " echo "${samplename} $outputdir/${b1}.paired.fq.gz $outputdir/${b2}.paired.fq.gz" > ${rootdir}/sample.information " >> $qsub1
	   `chmod g+r $qsub1 `
	   jobid=`qsub $qsub1`
	   `qhold -h u $jobid`
	   echo $jobid >> $TopOutputLogs/pbs.TRIM
           echo $jobid >> $TopOutputLogs/pbs.summary_dependencies_qc_trimming

	   echo `date`
	   `qrls -h u $jobid`
	  
    fi           
    done < $sampleinfo
	
	   jobids=$( cat $TopOutputLogs/pbs.summary_dependencies_qc_trimming | sed "s/\.[a-z]*//g" | tr "\n" ":" )
           qsub1=$TopOutputLogs/qsub.update.sampleinfo
	   cat $generic_qsub_header > $qsub1
           echo "#PBS -N update.sampleinfo" >> $qsub1
           echo "#PBS -e $TopOutputLogs/qsub.trim.samples.er" >> $qsub1
           echo "#PBS -o $TopOutputLogs/qsub.trim.samples.ou" >> $qsub1
	   echo "#PBS -W depend=afterok:$jobids" >> $qsub1
           echo "set -x" >> $qsub1
           echo "Updating the SAMPLEINFORMATION file with the trimmed reads" >> $qsub1
           echo "mv ${rootdir}/sample.information ${sampleinfo}" >> $qsub1
	   jobid=`qsub $qsub1`
	   echo $jobid >> $TopOutputLogs/pbs.summary_dependencies_qc_trimming
	   echo `date`
fi

MSG="Quality control for all samples finished successfully"
echo -e "program=$0 at line=$LINENO.\nReason=$MSG\n Logs for the run are available in $TopOutputLogs" | mail -s "[Task #${reportticket}]" "$redmine,$email"

echo `date`

set +x
echo -e "\n\n##################################################################################" >&2
echo -e "#############    DONE PROCESSING SAMPLE $SampleName. EXITING NOW.  ###############" >&2
echo -e "##################################################################################\n\n" >&2
set -x

