#!/bin/bash

#This scrips performs a parameter sweep for the bqsr stage of a variant calling pipeline. Run it as: bqsr_sweep.sh <runfile> <sample> <chr> <log.in> <log.ou> <qsub>
set -x
redmine=hpcbio-redmine@igb.illinois.edu

if [ $# != 6 ]
then
     MSG="Parameter mismatch.\nRerun like this: $0 <runfile> <sample> <chr> <log.in> <log.ou> <qsub>\n"     
     echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "ParameterSweep for the BQSR stage failed" "$redmine"
     exit 1;
fi

echo -e "\n\n#####################################################################################"        
echo -e "#############             BEGIN ANALYSIS PROCEDURE                    ###############"
echo -e "#####################################################################################\n\n"        

echo -e "\n\n#####################################################################################"        
echo -e "#############             DECLARING VARIABLES                         ###############"
echo -e "#####################################################################################\n\n"        

scriptfile=$0
runfile=$1
SampleName=$2
chr=$3
elog=$4
olog=$5
qsubfile=$6
LOGS="jobid:${PBS_JOBID}\nqsubfile=$qsubfile\nerrorlog=$elog\noutputlog=$olog"

if [ ! -s $runfile ] ; then
	    MSG="$runfile runfile not found"
	    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "Variant Calling Workflow failure message" "$redmine"
	    exit 1;
fi

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
echo -e "########### PREP WORK: split aligned.bam by chr and grab indels for chr            ###" 
echo -e "##################################################################################\n\n"

cd $RealignDir

$samtoolsdir/samtools view -bu -@ $thr -h $inputbam $chr > $dedupsortedbam
exitcode=$?

echo `date`
if [ $exitcode -ne 0 ]; then
	MSG="samtools command failed to split $inputbam by chr exitcode=$exitcode"
	echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
	exit $exitcode;
fi
if [ -s $dedupsortedbam ]; then
    set +x         
    echo -e "### the file was created. But we are not done.     #############"
    echo -e "### sometimes we may have a BAM file with NO alignmnets      ###"
    set -x
    numAlignments=$( $samtoolsdir/samtools view -c $dedupsortedbam )

    echo `date`
    if [ $numAlignments -eq 0 ]
    then
        echo -e "${SampleName}\tREALIGNMENT\tFAIL\tsamtools command did not produce alignments for $dedupsortedbam\n" >> $qcfile
        MSG="samtools command did not produce alignments for $dedupsortedbam"
        echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
        exit 1;
    else
        set +x
        echo -e "####### $dedupsortedbam seems to be in order ###########"
        set -x 
    fi
else
    MSG="samtools command did not produce a file $dedupsortedbam"
    echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

echo -e "\n### index  $dedupsortedbam                                 ###\n"  

$samtoolsdir/samtools index $dedupsortedbam
exitcode=$?
echo `date`
if [ $exitcode -ne 0 ]
then
         MSG="samtools index command failed to split $inputbam by chr exitcode=$exitcode"
         echo -e "program=$scriptfile stopped at line=$LINENO.\nReason=$MSG\n$LOGS"
         exit $exitcode;
fi
set +x
echo -e "\n\n### grab indels for this region. We need two variables for recalibration"
echo -e "###  because they are specified differently for each GATK tool                            ##############\n\n"
set -x

cd $indel_local

recalparmsindels=$( find ${PWD} -name "*${chr}.vcf" | sed "s/^/ --knownSites /g" | tr "\n" " " )
recalparmsdbsnp=" -knownSites $dbsnp_local "

if [ `expr ${#recalparmsindels}` -lt 1 ] 
then
	MSG="no indels were found for $chr in this folder $indel_local"
	echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
	exit 1;
fi

echo -e "\n\n##################################################################################"
echo -e "################# Running BQSR: Default parameters            ##########################" 
echo -e "\n\n##################################################################################"

if [ ! -d  $rootdir/sweepBQSR.$SampleName/default ]; then
	mkdir -p  $rootdir/sweepBQSR.$SampleName/default
fi

cd $rootdir/sweepBQSR.$SampleName/default

set +x
module load R ################################################ R is assumed availble in the user's environment, along with some libraries. More details are on: https://software.broadinstitute.org/gatk/guide/article?id=2899
set -x

echo Run using the default parameters: > "../bqsr.summary.txt"
echo "The parameters are:     ( ics maxCycle mcs bqsrBAQGOP ddq idq lqt mdq )">> "../bqsr.summary.txt"
echo "The default values are: ( 3     500     2      40     45   45  2  -1  )">> "../bqsr.summary.txt"
echo  "------------------------------------------------------------------------------------------------------------" >> "../bqsr.summary.txt"

START=$(date +%s)
$javadir/java -Xmx16g  -Djava.io.tmpdir=$tmpdir -jar $gatkdir/GenomeAnalysisTK.jar\
	-T BaseRecalibrator\
	-R $ref_local\
	-I ${RealignDir}/${dedupsortedbam}\
	$recalparmsindels \
       	$recalparmsdbsnp \
	--out $SampleName.$chr.recal_table.default.0 \
	-nct $thr

END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from BaseRecalibrator chr=$chr using: defaults 0 = $exitcode ; execution time = $DIFF" >> "../bqsr.summary.txt"
	
START=$(date +%s)
$javadir/java -Xmx16g  -Djava.io.tmpdir=$tmpdir -jar $gatkdir/GenomeAnalysisTK.jar\
       	-T BaseRecalibrator\
        -R $ref_local\
	-I ${RealignDir}/${dedupsortedbam}\
        $recalparmsindels \
	$recalparmsdbsnp \
        -BQSR $SampleName.$chr.recal_table.default.0\
        -o $SampleName.$chr.after_recal_table.default.0\
	-nct $thr
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from the rerun of BaseRecalibrator chr=$chr using: defaults 0 = $exitcode  ; execution time = $DIFF" >> "../bqsr.summary.txt"

START=$(date +%s)
$javadir/java -Xmx16g  -Djava.io.tmpdir=$tmpdir -jar $gatkdir/GenomeAnalysisTK.jar\
	-T AnalyzeCovariates\
	-R $ref_local\
	-before $SampleName.$chr.recal_table.default.0\
	-after $SampleName.$chr.after_recal_table.default.0\
	-plots $SampleName.$chr.recalibration_plots.default.0.pdf
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from AnalyzeCovariates chr=$chr using: defaults 0 = $exitcode  ; execution time = $DIFF" >> "../bqsr.summary.txt"
		
START=$(date +%s)
$javadir/java -Xmx16g  -Djava.io.tmpdir=$tmpdir -jar $gatkdir/GenomeAnalysisTK.jar\
	-T PrintReads\
	-R $ref_local\
	-I ${RealignDir}/${dedupsortedbam}\
	-BQSR $SampleName.$chr.recal_table.default.0\
	-o $SampleName.$chr.recal.default.0.bam\
	-nct $thr
END=$(date +%s)
DIFF=$(( $END - $START ))
exitcode=$?
echo "exit code from PrintReads chr=$chr using: defaults 0 = $exitcode  ; execution time = $DIFF" >> "../bqsr.summary.txt"
	
echo >> "../bqsr.summary.txt"
echo "####################################################################################################" >> ../bqsr.summary.txt


declare -a parameters=(ics maxCycle mcs bqsrBAQGOP ddq idq lqt mdq) 
declare -a min=(1 250 1 10 10 10 1 2) 
declare -a step=(2 150 2 10 10 10 2 2) 
declare -a max=(13 1000 13 70 70 70 12 12)


cd $rootdir/sweepBQSR.$SampleName
mkdir ${parameters[@]}

echo The parameters being tested and their ranges are given below: >> "bqsr.summary.txt"
echo paramters: ${parameters[@]}, >> "bqsr.summary.txt"
echo minimum  : ${min[@]} >> "bqsr.summary.txt"
echo maximum  : ${max[@]} >> "bqsr.summary.txt"
echo >> "bqsr.summary.txt"
echo  "------------------------------------------------------------------------------------------------------------" >> "bqsr.summary.txt"

pos=0

while [ $pos -lt ${#parameters[@]} ]; do
	par=${parameters[pos]}
	cd $rootdir/sweepBQSR.$SampleName/$par
	for i in $(seq ${min[pos]} ${step[pos]} ${max[pos]}); do
		START=$(date +%s)
		$javadir/java -Xmx16g  -Djava.io.tmpdir=$tmpdir -jar $gatkdir/GenomeAnalysisTK.jar \
		        -T BaseRecalibrator\
		        -R $ref_local\
			-I ${RealignDir}/${dedupsortedbam}\
                        $recalparmsindels \
                	$recalparmsdbsnp \
			-$par "$i"\
		        -o $SampleName.$chr.recal_table.$par.$i\
			-nct $thr
		END=$(date +%s)
		DIFF=$(( $END - $START ))
		exitcode=$?
		echo "exit code from BaseRecalibrator using: -$par $i = $exitcode ; execution time = $DIFF" >> "../bqsr.summary.txt"

		if [ -s "$SampleName.$chr.recal_table.$par.$i" ]; then
			echo 'successful parameter test chr=$chr $par=$i'
		else
			echo "BQSR chr=$chr failed using -$par $i">> "../bqsr.summary.txt"
			echo "Execution time is :$DIFF: seconds" >> "../bqsr.summary.txt"
			echo
			continue
		fi
			
		START=$(date +%s)
		$javadir/java -Xmx16g  -Djava.io.tmpdir=$tmpdir -jar $gatkdir/GenomeAnalysisTK.jar\
		        -T BaseRecalibrator\
		        -R $ref_local\
			-I ${RealignDir}/${dedupsortedbam}\
                        $recalparmsindels \
                	$recalparmsdbsnp \
		        -BQSR $SampleName.$chr.recal_table.$par.$i\
		        -o $SampleName.$chr.after_recal_table.$par.$i\
			-nct $thr
		END=$(date +%s)
		DIFF=$(( $END - $START ))
		exitcode=$?
	        echo "exit code from the rerun of  BaseRecalibrator with: -$par $i = $exitcode ; execution time = $DIFF" >> "../bqsr.summary.txt"
		
		START=$(date +%s)
		$javadir/java -Xmx16g  -Djava.io.tmpdir=$tmpdir -jar $gatkdir/GenomeAnalysisTK.jar\
		        -T AnalyzeCovariates\
		        -R $ref_local\
		        -before $SampleName.$chr.recal_table.$par.$i\
		        -after $SampleName.$chr.after_recal_table.$par.$i\
		        -plots $SampleName.$chr.recalibration_plots.$par.$i.pdf
		END=$(date +%s)
		DIFF=$(( $END - $START ))
		exitcode=$?
	        echo "exit code from AnalyzeCovariates using $par $i = $exitcode; execution time = $DIFF" >> "../bqsr.summary.txt"

		START=$(date +%s)
		$javadir/java -Xmx16g  -Djava.io.tmpdir=$tmpdir -jar $gatkdir/GenomeAnalysisTK.jar\
		        -T PrintReads\
		        -R $ref_local\
			-I ${RealignDir}/${dedupsortedbam}\
		        -BQSR $SampleName.$chr.recal_table.$par.$i\
		        -o $SampleName.$chr.recal.$par.$i.bam\
			-nct $thr
		END=$(date +%s)
		DIFF=$(( $END -$START + $DIFF ))
		exitcode=$?
		echo "exit code from PrintReads chr=$chr with: -$par $i = $exitcode  ; execution time = $DIFF" >> "../bqsr.summary.txt"
		echo  >> "../bqsr.summary.txt"
	done
let pos+=1
done

echo -e "\n\n########################################################################################"
echo -e "##############                 plotting now for this chr=$chr                            ##################"	
echo -e "########################################################################################\n\n"

Rscript $scriptdir/ParameterSweep/bqsrplotting.R $rootdir/sweepBQSR.$SampleName/ $thr $chr



echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
