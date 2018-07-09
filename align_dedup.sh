#!/bin/bash

#---------------------------------------------------------------------------------------------------------------------------------
#
# align_dedup_MuTect.sh <runfile> <SampleName> <tumor_read1> <tumor_read2> <normal_read1> <normal_read2> <log.in> <log.ou> <qsub>
#
#---------------------------------------------------------------------------------------------------------------------------------





#-----------------------------------------------------------------------------------------------
## COMMAND CHECK
#-----------------------------------------------------------------------------------------------

redmine=hpcbio-redmine@igb.illinois.edu
##redmine=mrweber2@illinois.edu

ERRLOG=error.log

set -x

if [[ $# != 8 ]]
then
        MSG="Parameter mismatch. Rerun as: $0 <runfile> <SampleName> <tumor_read1> <tumor_read2> <normal_read1> <normal_read2> <log.in> <log.ou> <qsub>"
        echo -e "program=$0 stopped at line=${LINENO}. Reason=${MSG}" >> ${ERRLOG}
        exit 1;
fi

#-----------------------------------------------------------------------------------------------





#-----------------------------------------------------------------------------------------------
## DECLARE VARIABLES
#-----------------------------------------------------------------------------------------------

echo `date`

scriptfile=$0
runfile=$1
SampleName=$2
RT1=$3
RT2=$4
RN1=$5
RN2=$6
log=$7
command=$8
LOGS="scriptfile=${scriptfile}\nlog=${log}\n"

if [[ ! -s ${runfile} ]]
then
	MSG="${runfile} runfile not found"
	echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" >> ${ERRLOG}
	exit 1;
fi

reportticket=$( cat ${runfile} | grep -w REPORTTICKET | cut -d '=' -f2 )
rootdir=$( cat ${runfile} | grep -w OUTPUTDIR | cut -d '=' -f2 )
deliverydir=$( cat ${runfile} | grep -w DELIVERYFOLDER | cut -d '=' -f2 ) 
tmpdir=$( cat ${runfile} | grep -w TMPDIR | cut -d '=' -f2 )
analysis=$( cat ${runfile} | grep -w ANALYSIS | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
thr=$( cat ${runfile} | grep -w PBSCORES | cut -d '=' -f2 )
indeldir=$( cat ${runfile} | grep -w INDELDIR | cut -d '=' -f2 )
alignertool=$( cat ${runfile} | grep -w ALIGNERTOOL | cut -d '=' -f2  )
bwamemdir=$( cat ${runfile} | grep -w BWAMEMDIR | cut -d '=' -f2  )
novocraftdir=$( cat ${runfile} | grep -w NOVOCRAFTDIR | cut -d '=' -f2  )
bwamem_parms=$( cat ${runfile} | grep -w BWAMEMPARAMS | cut -d '=' -f2 )
novoalign_parms=$( cat ${runfile} | grep -w NOVOALIGNPARAMS | cut -d '=' -f2 )
bwa_index=$( cat ${runfile} | grep -w BWAINDEX | cut -d '=' -f2 )
novoalign_index=$( cat ${runfile} | grep -w NOVOALIGNINDEX | cut -d '=' -f2 )
samblasterdir=$( cat ${runfile} | grep -w SAMBLASTERDIR | cut -d '=' -f2 )
samtoolsdir=$( cat ${runfile} | grep -w SAMDIR | cut -d '=' -f2 )
markduplicates=$( cat ${runfile} | grep -w MARKDUPLICATESTOOL | cut -d '=' -f2 | tr '[a-z]' '[A-Z]' )
picardir=$( cat ${runfile} | grep -w PICARDIR | cut -d '=' -f2 )
javadir=$( cat ${runfile} | grep -w JAVADIR | cut -d '=' -f2 )
sPL=$( cat ${runfile} | grep -w SAMPLEPL | cut -d '=' -f2 )
sCN=$( cat ${runfile} | grep -w SAMPLECN | cut -d '=' -f2 )
sLB=$( cat ${runfile} | grep -w SAMPLELB | cut -d '=' -f2 )
dup_cutoff=$( cat ${runfile} | grep -w  DUP_CUTOFF | cut -d '=' -f2 )
map_cutoff=$( cat ${runfile} | grep -w  MAP_CUTOFF | cut -d '=' -f2 )
outputdir=${rootdir}/${SampleName}

echo -e "\n\n#-----------------------------------------------------------------------------------------------\n#####   Setting variables for full workflow...\n#-----------------------------------------------------------------------------------------------\n"

SampleDir=${outputdir}
AlignDir=${outputdir}/align
RealignDir=${outputdir}/recal
VarcallDir=${outputdir}/variant
deliverydir=${rootdir}/${deliverydir}/${SampleName}
qctumorfile=${rootdir}/${deliverydir}/docs/QC_tumor_test_results.txt            # name of the txt file with all QC test results
qcnormalfile=${rootdir}/${deliverydir}/docs/QC_normal_test_results.txt            # name of the txt file with all QC test results
alignedtumorbam=${SampleName}.tumor.nodups.bam                              # name of the aligned bam
alignednormalbam=${SampleName}.nodups.bam                              # name of the aligned bam
alignedsortedtumorbam=${SampleName}.tumor.nodups.sorted.bam                 # name of the aligned-sorted bam
alignedsortednormalbam=${SampleName}.nodups.sorted.bam                 # name of the aligned-sorted bam
deduptumorbam=${SampleName}.tumor.wdups.bam                                 # name of the deduplicated bam
dedupnormalbam=${SampleName}.wdups.bam
dedupsortedtumorbam=${SampleName}.tumor.wdups.sorted.bam                    # name of the dedup-sorted bam (output of this module)
dedupsortednormalbam=${SampleName}.wdups.sorted.bam                    # name of the dedup-sorted bam (output of this module)

#-----------------------------------------------------------------------------------------------





#-----------------------------------------------------------------------------------------------
## SANITY CHECK
#-----------------------------------------------------------------------------------------------
 
if [[ ! -d ${tmpdir} ]]
then
	mkdir -p ${tmpdir}
fi

if [[ ! -d ${rootdir} ]]
then
	MSG="Invalid value specified for OUTPUTDIR=${rootdir} in the runfile."
	echo -e "program=$0 stopped at line=${LINENO}. Reason=${MSG}" >> ${ERRLOG}
	exit 1;
fi

if [[ ! -d ${outputdir} ]]
then
	MSG="${outputdir} outputdir not found"
	echo -e "program=$0 stopped at line=${LINENO}. Reason=${MSG}" >> ${ERRLOG}
	exit 1;
fi

if [[ ! -d ${deliverydir} ]]
then
	mkdir -p ${deliverydir}
fi

if [[ `expr ${#RT1}` -lt 1 ]]
then
	MSG="${RT1} tumor read one file not found"
	echo -e "Program $0 stopped at line=${LINENO}.\n\n${MSG}" >> ${ERRLOG}
	exit 1;

elif [[ ! -s ${RT1} ]]
then
	MSG="${RT1} tumor read one file not found"
	echo -e "Program $0 stopped at line=${LINENO}.\n\n${MSG}" >> ${ERRLOG}
	exit 1;                
fi

if [[ `expr ${#RT2}` -lt 1 ]]
then
	MSG="${RT2} tumor read two file not found"
	echo -e "Program $0 stopped at line=${LINENO}.\n\n${MSG}" >> ${ERRLOG}
	exit 1;

elif [[ ! -s ${RT2} ]]
then
	MSG="${RT2} tumor read two  file not found"
	echo -e "Program $0 stopped at line=${LINENO}.\n\n${MSG}" >> ${ERRLOG}
	exit 1;
fi

if [[ `expr ${#RN1}` -lt 1 ]]
then
	MSG="${RN1} normal read one file not found"
	echo -e "Program $0 stopped at line=${LINENO}.\n\n${MSG}" >> ${ERRLOG}
	exit 1;

elif [[ ! -s ${RN1} ]]
then
	MSG="${RN1} normal read one file not found"
	echo -e "Program $0 stopped at line=${LINENO}.\n\n${MSG}" >> ${ERRLOG}
	exit 1;
fi

if [[ `expr ${#RN2}` -lt 1 ]]
then
	MSG="${RN2} normal read two file not found"
	echo -e "Program $0 stopped at line=${LINENO}.\n\n${MSG}" >> ${ERRLOG}
	exit 1;

elif [[ ! -s ${RN2} ]]
then
	MSG="${RN2} normal read two  file not found"
	echo -e "Program $0 stopped at line=${LINENO}.\n\n${MSG}" >> ${ERRLOG}
	exit 1;
fi

if [[ `expr ${#SampleName}` -lt 1 ]]
then
	MSG="${SampleName} sample undefined variable"
	echo -e "Program $0 stopped at line=${LINENO}.\n\n${MSG}" >> ${ERRLOG}
	exit 1;

else

	sID=${SampleName}
	sPU=${SampleName}
	sSM=${SampleName}

fi

if [[ `expr ${#sLB}` -lt 1 -o `expr ${#sPL}` -lt 1 -o `expr ${#sCN}` -lt 1 ]]
then
	MSG="SAMPLELB=$sLB SAMPLEPL=$sPL SAMPLECN=$sCN at least one of these fields has invalid values. "
	echo -e "program=$0 stopped at line=${LINENO}.\nReason=${MSG}" >> ${ERRLOG}
	exit 1;
fi

RGparms=$( echo "ID=${sID}:LB=${sLB}:PL=${sPL}:PU=${sPU}:SM=${sSM}:CN=${sCN}" )
rgheader=$( echo -n -e "@RG\t" )$( echo -e "${RGparms}"  | tr ":" "\t" | tr "=" ":" )

if [[ `expr ${#markduplicates}` -lt 1 ]]
then
	markduplicates="NOVOSORT"
fi

#-----------------------------------------------------------------------------------------------





#-----------------------------------------------------------------------------------------------
## ALIGN-DEDUPLICATION STAGE
#-----------------------------------------------------------------------------------------------
       
echo -e "\n\n#-----------------------------------------------------------------------------------------------\n#####   ALIGN-DEDUPLICATION FOR SAMPLE ${SampleName}\n#-----------------------------------------------------------------------------------------------\n\n"          

echo `date` 

cd  ${AlignDir}

if [[ ${markduplicates} == "SAMBLASTER" ]]
then
	echo -e "\n\n#-----------------------------------------------------------------------------------------------\n#####   CASE1: dedup tool is ${markduplicates}\n#####   Step One: Alignment and Deduplication\n\n"
	

	$bwamemdir/bwa mem $bwamem_parms -t $thr -R "${rgheader}" $bwa_index ${RT1} ${RT2} | $samblaster | $samtoolsdir/samtools view -@ $thr -bSu -> $deduptumorbam 
	exitcode=$?
	echo `date`

	if [[ $exitcode -ne 0 ]]
	then
		MSG="ALIGNMENT-DEDUPLICATION failed with exitcode=$exitcode for tumor ${SampleName}"
		echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
		exit $exitcode;
	fi

	`chmod 660 ${deduptumorbam}*`

	$bwamemdir/bwa mem $bwamem_parms -t $thr -R "${rgheader}" $bwa_index ${RN1} ${RN2} | $samblaster | $samtoolsdir/samtools view -@ $thr -bSu -> $dedupnormalbam
	exitcode=$?
	echo `date`

	if [[ $exitcode -ne 0 ]]
	then
		MSG="ALIGNMENT-DEDUPLICATION failed with exitcode=$exitcode for normal ${SampleName}"
		echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
		exit $exitcode;
	fi

	`chmod 660 ${dedupnormalbam}*`

	echo -e "#####   Step Two: Making sure that a file was produced with alignments\n\n"

	if [[ -s ${AlignDir}/$deduptumorbam ]]
	then 
		echo -e "#####   The file was created. But we are not done.\n##### Sometimes we may have a BAM file with NO alignmnets"

		numAlignments=$( $samtoolsdir/samtools view -c ${AlignDir}/$deduptumorbam )
		echo `date`
		if [[ $numAlignments -eq 0 ]]
		then
			echo -e "${SampleName}\tALIGNMENT\tFAIL\tbwa mem command did not produce alignments for ${AlignDir}/$deduptumorbam\n" >> $qctumorfile
			MSG="bwa mem command did not produce alignments for ${AlignDir}/$deduptumorbam"
			echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
			exit 1;
		else
			echo -e "#####   ${AlignDir}/${deduptumorbam} seems to be in order\n\n"
		fi
	else 
		MSG="bwa mem command did not produce a file ${AlignDir}/$deduptumorbam"
        	echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
		exit 1;          
	fi

	if [[ -s ${AlignDir}/$dedupnormalbam ]]
        then
		echo -e "#####   The file was created. But we are not done.\n#####   Sometimes we may have a BAM file with NO alignments\n\n"
		numAlignments=$( $samtoolsdir/samtools view -c ${AlignDir}/$dedupnormalbam )
		echo `date`
		if [[ $numAlignments -eq 0 ]]
		then
			echo -e "${SampleName}\tALIGNMENT\tFAIL\tbwa mem command did not produce alignments for ${AlignDir}/$dedupnormalbam\n" >> $qcnormalfile
			MSG="bwa mem command did not produce alignments for ${AlignDir}/$dedupnormalbam"
			echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
			exit 1;
		else
			echo -e "#####   ${AlignDir}/$dedupnormalbam seems to be in order \n\n"
		fi
	else
		MSG="bwa mem command did not produce a file ${AlignDir}/$dedupnormalbam"
		echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
		exit 1;
	fi

	echo -e "#####   Step Three: Sort\n\n"

	$novocraftdir/novosort --index --tmpdir ${tmpdir} --threads $thr --compression 1 -o $dedupsortedtumorbam $deduptumorbam
	exitcode=$?

	`chmod 660 ${dedupsortedtumorbam}*`
	echo `date`

	if [[ $exitcode -ne 0 ]]
	then
		MSG="align-sorting step failed for tumor sample ${SampleName} exitcode=$exitcode."
		echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" >> ${ERRLOG}
		exit $exitcode;
	fi
	
	$novocraftdir/novosort --index --tmpdir ${tmpdir} --threads $thr --compression 1 -o $dedupsortednormalbam $dedupnormalbam
	exitcode=$?

	`chmod 660 ${dedupsortednormalbam}*`
	echo `date`
        
	if [[ $exitcode -ne 0 ]]
	then
		MSG="align-sorting step failed for normal sample ${SampleName} exitcode=$exitcode."
		echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" >> ${ERRLOG}
		exit $exitcode;
	fi

	echo -e "#####   Step Four: Making sure that a file was produced with alignments\n\n"

	if [[ -s ${AlignDir}/$dedupsortedtumorbam ]]
	then 
		echo -e "#####   The file was created. But we are not done.\n#####   Sometimes we may have a BAM file with NO alignments\n\n"
		numAlignments=$( $samtoolsdir/samtools view -c ${AlignDir}/$dedupsortedtumorbam ) 

		echo `date`
		if [[ $numAlignments -eq 0 ]]
		then
			echo -e "${SampleName}\tALIGNMENT\tFAIL\tnovosort command did not produce a file for ${AlignDir}/$dedupsortedtumorbam\n" >> $qctumorfile
			MSG="novosort command did not produce a file for ${AlignDir}/$dedupsortedtumorbam"
			echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
			exit 1;
		else
			echo -e "#####   ${AlignDir}/$deduptumorbam seems to be in order\n\n"
       
		fi
	else 
		MSG="novosort command did not produce a file ${AlignDir}/$dedupsortedtumorbam"
		echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
		exit 1;          
	fi       	

	if [[ -s ${AlignDir}/$dedupsortednormalbam ]]
	then
		echo -e "#####   The file was created. But we are not done.\n#####   Sometimes we may have a BAM file with NO alignmnets\n\n" 
		numAlignments=$( $samtoolsdir/samtools view -c ${AlignDir}/$dedupsortednormalbam )

		echo `date`
		if [[ $numAlignments -eq 0 ]]
		then
			echo -e "${SampleName}\tALIGNMENT\tFAIL\tnovosort command did not produce a file for ${AlignDir}/$dedupsortednormalbam\n" >> $qcnormalfile
			MSG="novosort command did not produce a file for ${AlignDir}/$dedupsortednormalbam"
			echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
			exit 1;
		else
			echo -e "#####   ${AlignDir}/$dedupnormalbam seems to be in order\n\n"

		fi
	else
		MSG="novosort command did not produce a file ${AlignDir}/$dedupsortednormalbam"
		echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
		exit 1;
	fi

elif  [[ ${markduplicates} == "NOVOSORT" ]]
then
	echo -e "\n\n#####   CASE2: dedup tool is NOVOSORT. one cmd for align and one for dedup-sort\n\n#####    Step One: Alignment\n\n" 

	if [[ $alignertool== "BWA" ]]
	then

		$bwamemdir/bwa mem $bwamem_parms -t $thr -R "${rgheader}" $bwa_index ${RT1} ${RT2} | $samtoolsdir/samtools view -@ $thr -bSu -> $alignedtumorbam 
		exitcode=$?
        	`chmod 660 ${alignedtumorbam}*`
		echo `date
`
		if [[ $exitcode -ne 0 ]]
		then
	 		MSG="alignment step  failed for tumor sample ${SampleName} exitcode=$exitcode."
	 		echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" >> ${ERRLOG}
			exit $exitcode;
		fi

	elif [[ $alignertool == "NOVOALIGN" ]]
	then
		$novocraftdir/novoalign $novoalign_parms  -c $thr -d ${novoalign_index} -f ${RT1} ${RT2} -o SAM | $samtoolsdir/samtools view -@ $thr -bS - > $alignedtumorbam
		exitcode=$?
		`chmod 660 ${alignedtumorbam}*`
		
		echo `date`
		
		if [[ $exitcode -ne 0 ]]
		then
			MSG="alignment step  failed for tumor sample ${SampleName} exitcode=$exitcode."
			echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" >> ${ERRLOG}
			exit $exitcode;
		fi

		$bwamemdir/bwa mem $bwamem_parms -t $thr -R "${rgheader}" $bwa_index ${RN1} ${RN2} | $samtoolsdir/samtools view -@ $thr -bSu -> $alignednormalbam
		exitcode=$?
		`chmod 660 ${alignednormalbam}*`
		
		echo `date`
           
		if [[ $exitcode -ne 0 ]]
		then
			MSG="alignment step  failed for normal sample ${SampleName} exitcode=$exitcode."
			echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" >> ${ERRLOG}
			exit $exitcode;
		fi
	elif [[ $alignertool == "NOVOALIGN" ]]
	then
		$novocraftdir/novoalign $novoalign_parms  -c $thr -d ${novoalign_index} -f ${RN1} ${RN2} -o SAM | $samtoolsdir/samtools view -@ $thr -bS - > $alignednormalbam
		exitcode=$?
		`chmod 660 ${alignednormalbam}*`
           
		echo `date`
           
		if [[ $exitcode -ne 0 ]]
		then
			MSG="alignment step  failed for normal sample ${SampleName} exitcode=$exitcode."
			echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" >> ${ERRLOG}
			exit $exitcode;
		fi
	fi   
	 
	echo -e "#####   Step Two: Making sure that a file was produced with alignments\n\n"

	if [[ -s ${AlignDir}/$alignedtumorbam ]]
	then
		echo -e "#####   The file was created, but we are not done.\n#####   Sometimes we may have a BAM file with NO alignments\n\n" 		

		numAlignments=$( $samtoolsdir/samtools view -c ${AlignDir}/$alignedtumorbam ) 

		echo `date`
		
		if [[ $numAlignments -eq 0 ]]
		then
			echo -e "${SampleName}\tALIGNMENT\tFAIL\taligner command did not produce alignments for ${AlignDir}/$alignedtumorbam\n" >> $qctumorfile
			MSG="aligner command did not produce alignments for ${AlignDir}/$alignedtumorbam"
			echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
			exit 1;
		else
			echo -e "#####   ${AlignDir}/$alignedtumorbam seems to be in order\n\n" 
		fi
	else 
		MSG="aligner command did not produce a file ${AlignDir}/$alignedtumorbam"
		echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
		exit 1;          
	fi 

	if [[ -s ${AlignDir}/$alignednormalbam ]]
	then
		echo -e "#####   The file was created, but we are not done.\n#####   Sometimes we may have a BAM file with NO alignments\n\n"              

		numAlignments=$( $samtoolsdir/samtools view -c ${AlignDir}/$alignednormalbam )

		echo `date`
		
		if [[ $numAlignments -eq 0 ]]
		then
			echo -e "${SampleName}\tALIGNMENT\tFAIL\taligner command did not produce alignments for ${AlignDir}/$alignednormalbam\n" >> $qcnormalfile
			MSG="aligner command did not produce alignments for ${AlignDir}/$alignednormalbam"
			echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
			exit 1;
		else
			echo -e "#####   ${AlignDir}/$alignednormalbam seems to be in order\n\n"
		fi
	else
		MSG="aligner command did not produce a file ${AlignDir}/$alignednormalbam"
		echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
		exit 1;
	fi

	echo -e "#####   Step Three: Sort + Dedup + Indexing\n\n"

	$novocraftdir/novosort -r "${rgheader}" --markDuplicates  -t ${tmpdir} -c $thr -i -o $dedupsortedtumorbam $alignedtumorbam
	exitcode=$?
	`chmod 660 ${dedupsortedtumorbam}*`
	
	echo `date`
	
	if [[ $exitcode -ne 0 ]]
	then
		MSG="alignment step  failed during sorting-deduplication for tumor sample ${SampleName} exitcode=$exitcode."
		echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" >> ${ERRLOG}
		exit $exitcode;
	fi
 
	$novocraftdir/novosort -r "${rgheader}" --markDuplicates  -t ${tmpdir} -c $thr -i -o $dedupsortednormalbam $alignednormalbam
	exitcode=$?
	`chmod 660 ${dedupsortednormalbam}*`
        
	echo `date`

	if [[ $exitcode -ne 0 ]]
	then
		MSG="alignment step  failed during sorting-deduplication for normal sample ${SampleName} exitcode=$exitcode."
		echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" >> ${ERRLOG}
		exit $exitcode;
	fi
	
	echo -e "#####   Step Four: Making sure that a file was produced with alignments\n\n" >&2
	
	if [[ -s ${AlignDir}/$dedupsortedtumorbam ]]
	then
		echo -e "#####   The file was created, but we are not done\n#####   Sometimes we may have a BAM file with NO alignments\n\n" 
	    
		numAlignments=$( $samtoolsdir/samtools view -c ${AlignDir}/$dedupsortedtumorbam ) 

		echo `date`
		if [[ $numAlignments -eq 0 ]]
		then
			echo -e "${SampleName}\tALIGNMENT\tFAIL\tnovosort command did not produce a file for ${AlignDir}/$dedupsortedtumorbam\n" >> $qctumorfile	    
			MSG="novosort command did not produce a file for ${AlignDir}/$dedupsortedtumorbam"
			echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
			exit 1;
		else
			echo -e "####### ${AlignDir}/$deduptumorbam seems to be in order\n\n"
		fi
	else
		MSG="novosort command did not produce a file ${AlignDir}/$dedupsortedtumorbam"
		echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
		exit 1;          
	fi   
	
	if [[ -s ${AlignDir}/$dedupsortednormalbam ]]
        then
		echo -e "#####   The file was created, but we are not done.\n\n#####   Sometimes we may have a BAM file with NO alignments\n\n" 

		numAlignments=$( $samtoolsdir/samtools view -c ${AlignDir}/$dedupsortednormalbam )

		echo `date`
		if [[ $numAlignments -eq 0 ]]
		then
			echo -e "${SampleName}\tALIGNMENT\tFAIL\tnovosort command did not produce a file for ${AlignDir}/$dedupsortednormalbam\n" >> $qcnormalfile
			MSG="novosort command did not produce a file for ${AlignDir}/$dedupsortednormalbam"
			echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
			exit 1;
		else
			echo -e "#####   ${AlignDir}/$dedupnormalbam seems to be in order\n\n"
		fi
	else
		MSG="novosort command did not produce a file ${AlignDir}/$dedupsortednormalbam"
		echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
		exit 1;
	fi

#-----------------------------------------------------------------------------------------------






	echo -e "\n\n##################################################################################"
	echo -e "#############      END NOVOSRT  BLOCK                                 ############"
	echo -e "##################################################################################\n\n"             
	

elif  [ ${markduplicates} == "PICARD" ]
then
	set +x
	echo -e "\n\n########################################################################################" >&2
	echo -e "####  CASE2: dedup tool is PICARD. one cmd for align, one for sort, one for dedup   ####" >&2 
	echo -e "########################################################################################\n\n" >&2

	echo -e "\n\n###############################################################################" >&2	     
	echo -e "#############  step one: alignment                                 ############" >&2
	echo -e "###############################################################################\n\n" >&2
	set -x 	

	$bwamemdir/bwa mem $bwamem_parms -t $thr -R "${rgheader}" $bwa_index ${RT1} ${RT2} | $samtoolsdir/samtools view -@ $thr -bSu -> $alignedtumorbam 
	exitcode=$?
        `chmod 660 ${alignedtumorbam}*`
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="alignment step  failed for tumor sample ${SampleName} exitcode=$exitcode."
	    echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" | mail -s "[Task #${reportticket}]" "${redmine},$email"
	    exit $exitcode;
	fi   

        $bwamemdir/bwa mem $bwamem_parms -t $thr -R "${rgheader}" $bwa_index ${RN1} ${RN2} | $samtoolsdir/samtools view -@ $thr -bSu -> $alignednormalbam
        exitcode=$?
        `chmod 660 ${alignednormalbam}*`
        echo `date`
        if [ $exitcode -ne 0 ]
        then
            MSG="alignment step  failed for normal sample ${SampleName} exitcode=$exitcode."
            echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" | mail -s "[Task #${reportticket}]" "${redmine},$email"
            exit $exitcode;
        fi


	set +x
	echo -e "\n\n######################################################################################" >&2	     
	echo -e "#############  step two:making sure that a file was produced with alignments     #####" >&2
	echo -e "######################################################################################\n\n" >&2
	set -x

	if [ -s ${AlignDir}/$alignedtumorbam ]
	then 
	    set +x		           
	    echo -e "### the file was created. But we are not done.     #############" >&2
	    echo -e "### sometimes we may have a BAM file with NO alignmnets      ###" >&2
	    set -x

	    numAlignments=$( $samtoolsdir/samtools view -c ${AlignDir}/$alignedtumorbam ) 

	    echo `date`
	    if [ $numAlignments -eq 0 ]
	    then
	        echo -e "${SampleName}\tALIGNMENT\tFAIL\tbwa mem command did not produce alignments for ${AlignDir}/$alignedtumorbam\n" >> $qctumorfile	    
		MSG="bwa mem command did not produce alignments for ${AlignDir}/$alignedtumorbam"
                echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
		exit 1;
	    else
		set +x
		echo -e "####### ${AlignDir}/$alignedtumorbam seems to be in order ###########" >&2
		set -x 
	    fi
	else 
	    MSG="bwa mem command did not produce a file ${AlignDir}/$alignedtumorbam"
            echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
	    exit 1;          
	fi      

        if [ -s ${AlignDir}/$alignednormalbam ]
        then
            set +x                         
            echo -e "### the file was created. But we are not done.     #############" >&2
            echo -e "### sometimes we may have a BAM file with NO alignmnets      ###" >&2
            set -x

            numAlignments=$( $samtoolsdir/samtools view -c ${AlignDir}/$alignednormalbam )

            echo `date`
            if [ $numAlignments -eq 0 ]
            then
                echo -e "${SampleName}\tALIGNMENT\tFAIL\tbwa mem command did not produce alignments for ${AlignDir}/$alignednormalbam\n" >> $qcnormalfile
                MSG="bwa mem command did not produce alignments for ${AlignDir}/$alignednormalbam"
                echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
                exit 1;
            else
                set +x
                echo -e "####### ${AlignDir}/$alignednormalbam seems to be in order ###########" >&2
                set -x 
            fi
        else
            MSG="bwa mem command did not produce a file ${AlignDir}/$alignednormalbam"
            echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
            exit 1;
        fi


	set +x 
	echo -e "\n\n###################################################################################" >&2
	echo -e "#############  step three: sort                                        ############" >&2
	echo -e "###################################################################################\n\n" >&2
	set -x

	$novocraftdir/novosort -t ${tmpdir} -c ${thr} -i -o $alignedsortedtumorbam $alignedtumorbam
	exitcode=$?
        `chmod 660 ${alignedsortedtumorbam}*`
	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="alignment step  failed during sorting for tumor sample ${SampleName} exitcode=$exitcode."
	    echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" | mail -s "[Task #${reportticket}]" "${redmine},$email"
            exit $exitcode
	fi 

        $novocraftdir/novosort -t ${tmpdir} -c ${thr} -i -o $alignedsortednormalbam $alignednormalbam
        exitcode=$?
        `chmod 660 ${alignedsortednormalbam}*`
        echo `date`
        if [ $exitcode -ne 0 ]
        then
            MSG="alignment step  failed during sorting for normal sample ${SampleName} exitcode=$exitcode."
            echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" | mail -s "[Task #${reportticket}]" "${redmine},$email"
            exit $exitcode
        fi

	set +x
	echo -e "\n\n###################################################################################" >&2
	echo -e "#############  step four: dedup                                        ############" >&2
	echo -e "###################################################################################\n\n" >&2
	set -x 


        $javadir/java -Xmx8g -Djava.io.tmpdir=${tmpdir} -jar $picardir/picard.jar  MarkDuplicates \
           INPUT=$alignedsortedtumorbam OUTPUT=$dedupsortedtumorbam TMP_DIR=${tmpdir} \
           ASSUME_SORTED=true MAX_RECORDS_IN_RAM=null CREATE_INDEX=true \
           METRICS_FILE=${SampleName}.picard.metrics \
           VALIDATION_STRINGENCY=SILENT

	exitcode=$?
        `find . -type d | xargs chmod -R 770`
        `find . -type f | xargs chmod -R 660`

	echo `date`
	if [ $exitcode -ne 0 ]
	then
	    MSG="alignment step  failed during deduplication for tumor sample ${SampleName} exitcode=$exitcode."
	    echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" | mail -s "[Task #${reportticket}]" "${redmine},$email"
            exit $exitcode
	fi 

        $javadir/java -Xmx8g -Djava.io.tmpdir=${tmpdir} -jar $picardir/picard.jar  MarkDuplicates \
           INPUT=$alignedsortednormalbam OUTPUT=$dedupsortednormalbam TMP_DIR=${tmpdir} \
           ASSUME_SORTED=true MAX_RECORDS_IN_RAM=null CREATE_INDEX=true \
           METRICS_FILE=${SampleName}.picard.metrics \
           VALIDATION_STRINGENCY=SILENT

        exitcode=$?
        `find . -type d | xargs chmod -R 770`
        `find . -type f | xargs chmod -R 660`

        echo `date`
        if [ $exitcode -ne 0 ]
        then
            MSG="alignment step  failed during deduplication for normal sample ${SampleName} exitcode=$exitcode."
            echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" | mail -s "[Task #${reportticket}]" "${redmine},$email"
            exit $exitcode
        fi


	set +x	
	echo -e "\n\n######################################################################################" >&2	     
	echo -e "#############  step five: making sure that a file was produced with alignments #######" >&2
	echo -e "######################################################################################\n\n" >&2
	set -x
	
	if [ -s ${AlignDir}/$dedupsortedtumorbam ]
	then
	    set +x			
	    echo -e "### the file was created. But we are not done.     #############" >&2
	    echo -e "### sometimes we may have a BAM file with NO alignmnets      ###" >&2
	    set -x 
	    
	    numAlignments=$( $samtoolsdir/samtools view -c ${AlignDir}/$dedupsortedtumorbam ) 

	    echo `date`
	    if [ $numAlignments -eq 0 ]
	    then
	        echo -e "${SampleName}\tALIGNMENT\tFAIL\tpicard command did not produce a file for ${AlignDir}/$dedupsortedtumorbam\n" >> $qctumorfile	    
		MSG="novosort command did not produce a file for ${AlignDir}/$dedupsortedtumorbam"
                echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
		exit 1;
	    else
		set +x
		echo -e "####### ${AlignDir}/$dedupsortedtumorbam seems to be in order ###########" >&2
		set -x
	    fi
	else 
	    MSG="picard command did not produce a file ${AlignDir}/$dedupsortedtumorbam"
            echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
	    exit 1;          
	fi   

        if [ -s ${AlignDir}/$dedupsortednormalbam ]
        then
            set +x                      
            echo -e "### the file was created. But we are not done.     #############" >&2
            echo -e "### sometimes we may have a BAM file with NO alignmnets      ###" >&2
            set -x 

            numAlignments=$( $samtoolsdir/samtools view -c ${AlignDir}/$dedupsortednormalbam )

            echo `date`
            if [ $numAlignments -eq 0 ]
            then
                echo -e "${SampleName}\tALIGNMENT\tFAIL\tpicard command did not produce a file for ${AlignDir}/$dedupsortednormalbam\n" >> $qcnormalfile
                MSG="novosort command did not produce a file for ${AlignDir}/$dedupsortednormalbam"
                echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
                exit 1;
            else
                set +x
                echo -e "####### ${AlignDir}/$dedupsortednormalbam seems to be in order ###########" >&2
                set -x
            fi
        else
            MSG="picard command did not produce a file ${AlignDir}/$dedupsortednormalbam"
            echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
            exit 1;
        fi


        set +x 
	echo -e "\n\n#################################################################################" >&2
	echo -e "#############      END PICARD  BLOCK                                 ############" >&2
	echo -e "#################################################################################\n\n" >&2
	set -x             

	
else
	MSG="unrecognized deduplication tool ${markduplicates}"
        echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.FAILURE
	exit 1;        

fi


set +x 
echo -e "\n\n##################################################################################" >&2
echo -e "#############     END ALIGNMENT-DEDUPLICATION BLOCK                   ############" >&2
echo -e "##################################################################################\n\n" >&2

echo `date`

echo -e "\n\n##################################################################################" >&2
echo -e "##################################################################################" >&2          
echo -e "##################################################################################" >&2        
echo -e "########   ALIGNMENT QC TEST   FOR SAMPLE ${SampleName}                             " >&2
echo -e "########   QC rule1: duplication cutoff <= $dup_cutoff                            " >&2
echo -e "########   QC rule2: mapped_reads cutoff >= $map_cutoff                           " >&2
echo -e "##################################################################################" >&2
echo -e "##################################################################################" >&2 
echo -e "##################################################################################\n\n" >&2
        
     

echo -e "\n\n##################################################################################" >&2	     
echo -e "#############  step one: generating the relevant file with flagstat   ############" >&2
echo -e "##################################################################################\n\n" >&2
set -x

flagstatstumor=${dedupsortedtumorbam}.flagstats

echo `date`             
$samtoolsdir/samtools flagstat $dedupsortedtumorbam > $flagstatstumor
echo `date`
`find . -type d | xargs chmod -R 770`
`find . -type f | xargs chmod -R 660`

flagstatsnormal=${dedupsortednormalbam}.flagstats

echo `date`             
$samtoolsdir/samtools flagstat $dedupsortednormalbam > $flagstatsnormal
echo `date`
`find . -type d | xargs chmod -R 770`
`find . -type f | xargs chmod -R 660`


set +x
echo -e "\n\n##################################################################################" >&2	     
echo -e "#############  step two: sanity check                                 ############" >&2
echo -e "##################################################################################\n\n" >&2
set -x  

if [ ! -s $flagstatstumor ]
then
	 echo -e "${SampleName}\tQCT/EST\tFAIL\tsamtools/samtools flagstat command produced an empty file $flagstatstumor\n" >> $qctumorfile
	 MSG="samtools flagstat command produced an empty file  $flagstatstumor"
	 echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" | mail -s "[Task #${reportticket}]" "${redmine},$email"
	 exit $exitcode;
fi

if [ ! -s $flagstatsnormal ]
then
         echo -e "${SampleName}\tQCT/EST\tFAIL\tsamtools/samtools flagstat command produced an empty file $flagstatsnormal\n" >> $qcnormalfile
         MSG="samtools flagstat command produced an empty file  $flagstatsnormal"
         echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" | mail -s "[Task #${reportticket}]" "${redmine},$email"
         exit $exitcode;
fi


set +x
echo -e "\n\n######################################################################################" >&2	     
echo -e "#############  step three: parsing the file and grabbing stats for the QC test     ###" >&2
echo -e "######################################################################################\n\n" >&2
set -x           


tot_tumormapped=$( cat $flagstatstumor | grep "mapped (" | cut -d ' ' -f1 )
tot_tumorreads=$( cat $flagstatstumor | grep "in total" | cut -d ' ' -f1 )
tot_tumordups=$( cat $flagstatstumor | grep "duplicates" | cut -d ' ' -f1 )

#now testing if these variables are numeric and have numbers

if [ $tot_tumordups -eq $tot_tumordups 2>/dev/null ]
then
	echo -e "ok val"
else
	MSG="$flagstatstumor samtools flagstat file parsed incorrectly"
	echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" | mail -s "[Task #${reportticket}]" "${redmine},$email"
	exit $exitcode;
fi
if [ $tot_tumorreads -eq $tot_tumorreads 2>/dev/null ]
then
	echo -e "ok val"
else
	MSG="$flagstatstumor samtools flagstat file parsed incorrectly"
	echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" | mail -s "[Task #${reportticket}]" "${redmine},$email"
	exit $exitcode;
fi

if [ $tot_tumormapped -eq $tot_tumormapped 2>/dev/null ]
then
	echo -e "ok val"
else
	MSG="$flagstatstumor samtools flagstat file parsed incorrectly"
	echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" | mail -s "[Task #${reportticket}]" "${redmine},$email"
	exit $exitcode;
fi


tot_normalmapped=$( cat $flagstatsnormal | grep "mapped (" | cut -d ' ' -f1 )
tot_normalreads=$( cat $flagstatsnormal | grep "in total" | cut -d ' ' -f1 )
tot_normaldups=$( cat $flagstatsnormal | grep "duplicates" | cut -d ' ' -f1 )

#now testing if these variables are numeric and have numbers

if [ $tot_normaldups -eq $tot_normaldups 2>/dev/null ]
then
        echo -e "ok val"
else
        MSG="$flagstatsnormal samtools flagstat file parsed incorrectly"
        echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" | mail -s "[Task #${reportticket}]" "${redmine},$email"
        exit $exitcode;
fi
if [ $tot_normalreads -eq $tot_normalreads 2>/dev/null ]
then
        echo -e "ok val"
else
        MSG="$flagstatsnormal samtools flagstat file parsed incorrectly"
        echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" | mail -s "[Task #${reportticket}]" "${redmine},$email"
        exit $exitcode;
fi

if [ $tot_normalmapped -eq $tot_normalmapped 2>/dev/null ]
then
        echo -e "ok val"
else
        MSG="$flagstatsnormal samtools flagstat file parsed incorrectly"
        echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" | mail -s "[Task #${reportticket}]" "${redmine},$email"
        exit $exitcode;
fi


set +x           
echo -e "\n\n##################################################################################" >&2	     
echo -e "#############  step four: calculating stats according to QC rules              ###" >&2
echo -e "##################################################################################\n\n" >&2
set -x        


perc_tumordup=$(( tot_tumordups * 100 / tot_tumorreads ))
perc_tumormapped=$(( tot_tumormapped * 100 / tot_tumorreads ))

#now testing if these variables have numbers

if [ $perc_tumordup -eq $perc_tumordup 2>/dev/null ]
then
	echo -e "ok val"
else
	MSG="$flagstatstumor samtools flagstat file parsed incorrectly"
	echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" | mail -s "[Task #${reportticket}]" "${redmine},$email"
	exit $exitcode;
fi

if [ $perc_tumormapped -eq $perc_tumormapped 2>/dev/null ]
then
	echo -e "ok val"
else
	MSG="$flagstatstumor samtools flagstat file parsed incorrectly"
	echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" | mail -s "[Task #${reportticket}]" "${redmine},$email"
	exit $exitcode;
fi

perc_normaldup=$(( tot_normaldups * 100 / tot_normalreads ))
perc_normalmapped=$(( tot_normalmapped * 100 / tot_normalreads ))

#now testing if these variables have numbers

if [ $perc_normaldup -eq $perc_normaldup 2>/dev/null ]
then
        echo -e "ok val"
else
        MSG="$flagstatsnormal samtools flagstat file parsed incorrectly"
        echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" | mail -s "[Task #${reportticket}]" "${redmine},$email"
        exit $exitcode;
fi

if [ $perc_normalmapped -eq $perc_normalmapped 2>/dev/null ]
then
        echo -e "ok val"
else
        MSG="$flagstatsnormal samtools flagstat file parsed incorrectly"
        echo -e "program=${scriptfile} stopped at line=${LINENO}.\nReason=${MSG}\n${LOGS}" | mail -s "[Task #${reportticket}]" "${redmine},$email"
        exit $exitcode;
fi

set +x
echo -e "\n\n##################################################################################" >&2	     
echo -e "#############  step five: applying the  QC rules                               ###" >&2
echo -e "##################################################################################\n\n" >&2
set -x         

if [ $perc_normaldup -lt $dup_cutoff ]
then
	echo -e "$#####  sample passed first filter percent_duplicates with value $perc_normaldup, maximum cutoff is $dup_cutoff"
	
	if [ $perc_normalmapped -gt $map_cutoff ]
	then
	        echo -e "##### $sample passed second filter map_cutoff with value $perc_normalmapped, minimum cutoff is $map_cutoff"	
	        echo -e "${SampleName}\tQCTEST\tPASS\trule1 ok: percent_duplication=$perc_normaldup <? duplication_cutoff=$dup_cutoff\trule2 ok: percent_mapped=$perc_normalmapped >? mapping_cutoff=$map_cutoff" >> $qcnormalfile
	else
	        echo -e "##### $sample DID NOT pass second filter map_cutoff with value $perc_normalmapped, minimum cutoff is $map_cutoff"	
	        echo -e "${SampleName}\tQCTEST\tFAIL\trule1 ok: percent_duplication=$perc_normaldup <? duplication_cutoff=$dup_cutoff\trule2 notok: percent_mapped=$perc_normalmapped >? mapping_cutoff=$map_cutoff" >> $qcnormalfile	          
	fi
else
	echo -e "$#####  sample DID NOT pass first filter percent_duplicates with value $perc_normaldup, maximum cutoff is $dup_cutoff"
	echo -e "${SampleName}\tQCTEST\tFAIL\trule1 not ok: percent_duplication=$perc_normaldup <? duplication_cutoff=$dup_cutoff\trule2 not evaluated: percent_mapped=$perc_normalmapped >? mapping_cutoff=$map_cutoff" >> $qcnormalfile
fi

if [ $perc_tumordup -lt $dup_cutoff ]
then
        echo -e "$#####  sample passed first filter percent_duplicates with value $perc_tumordup, maximum cutoff is $dup_cutoff"

        if [ $perc_tumormapped -gt $map_cutoff ]
        then
                echo -e "##### $sample passed second filter map_cutoff with value $perc_tumormapped, minimum cutoff is $map_cutoff"  
                echo -e "${SampleName}\tQCTEST\tPASS\trule1 ok: percent_duplication=$perc_tumordup <? duplication_cutoff=$dup_cutoff\trule2 ok: percent_mapped=$perc_tumormapped >? mapping_cutoff=$map_cutoff" >> $qctumorfile
        else
                echo -e "##### $sample DID NOT pass second filter map_cutoff with value $perc_tumormapped, minimum cutoff is $map_cutoff"    
                echo -e "${SampleName}\tQCTEST\tFAIL\trule1 ok: percent_duplication=$perc_tumordup <? duplication_cutoff=$dup_cutoff\trule2 notok: percent_mapped=$perc_tumormapped >? mapping_cutoff=$map_cutoff" >> $qctumorfile
        fi
else
        echo -e "$#####  sample DID NOT pass first filter percent_duplicates with value $perc_tumordup, maximum cutoff is $dup_cutoff"
        echo -e "${SampleName}\tQCTEST\tFAIL\trule1 not ok: percent_duplication=$perc_tumordup <? duplication_cutoff=$dup_cutoff\trule2 not evaluated: percent_mapped=$perc_tumormapped >? mapping_cutoff=$map_cutoff" >> $qctumorfile
fi

set +x
echo -e "\n\n##################################################################################" >&2
echo -e "#############       END QC TEST                                       ############" >&2        
echo -e "##################################################################################\n\n" >&2

echo `date`

echo -e "\n\n##################################################################################"  >&2
echo -e "##################################################################################" >&2
echo -e "##################################################################################" >&2        
echo -e "#############   WRAP UP                                               ############" >&2       
echo -e "##################################################################################" >&2
echo -e "##################################################################################" >&2 
echo -e "##################################################################################\n\n" >&2	

echo `date`
set -x

### perhaps this bam file is not necessary in the delivery folder           
### cp ${AlignDir}/${SampleName}.wdups.sorted.bam          ${deliverydir}   


MSG="ALIGNMENT-DEDUPLICATION for ${SampleName} finished successfully"
echo -e "${MSG}" >> ${rootdir}/logs/mail.${analysis}.SUCCESS 

echo `date`

set +x
echo -e "\n\n#################################################################################################" >&2
echo -e "#############    DONE WITH ALIGNMENT-DEDUPLICATION ON SAMPLE ${SampleName}. EXITING NOW.  " >&2
echo -e "#################################################################################################\n\n" >&2
set -x
