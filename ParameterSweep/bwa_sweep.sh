#!/bin/bash

## This script sweeps the parameters for the aligner BWA MEM, and allows comaprison between the quality of mapping in each case. It should be called as: bwa_sweep.sh <runfile> <SampleName> <read1> <read2> <log.in> <log.ou> 

set -x 

scriptfile=$0
runfile=$1
SampleName=$2
read1=$3
read2=$4
redmine=hpcbio-redmine@igb.illinois.edu

if [ $# != 6 ]
then
	MSG="Parameter mismatch.\nRerun like this: $0 <runfile> <SampleName> <read1> <read2> <log.in> <log.ou> <qsub>\n"
	echo -e "program=$0 stopped at line=$LINENO. Reason=$MSG" | mail -s "ParameterSweep for the Alignment stage (bwa mem) failed" "$redmine"
	exit 1;
fi

echo -e "\n\n########################################################################################"
echo -e "#############                Pipeline starts here!              ###############"
echo -e "########################################################################################\n\n"

set -x
echo `date`
refdir=$( cat $runfile | grep -w REFGENOMEDIR | cut -d '=' -f2 )
bwa_index=$( cat $runfile | grep -w BWAINDEX | cut -d '=' -f2 )
outputdir=$( cat $runfile | grep -w OUTPUTDIR | cut -d '=' -f2 )
bwamemdir=$( cat $runfile | grep -w BWAMEMDIR | cut -d '=' -f2 )
samtoolsdir=$( cat $runfile | grep -w SAMDIR | cut -d '=' -f2 )
thr=$( cat $runfile | grep -w PBSCORES | cut -d '=' -f2 )
sPL=$( cat $runfile | grep -w SAMPLEPL | cut -d '=' -f2 )
sCN=$( cat $runfile | grep -w SAMPLECN | cut -d '=' -f2 )
sLB=$( cat $runfile | grep -w SAMPLELB | cut -d '=' -f2 )


results=$outputdir/sweepBWA.$SampleName

## Read group info:
if [ `expr ${#SampleName}` -lt 1 ]
then
    MSG="$SampleName sample undefined variable"
    echo -e "Program $0 stopped at line=$LINENO.\n\n$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1
else
    sID=$SampleName
    sPU=$SampleName
    sSM=$SampleName
fi
if [ `expr ${#sLB}` -lt 1 -o `expr ${#sPL}` -lt 1 -o `expr ${#sCN}` -lt 1 ]
then
    MSG="SAMPLELB=$sLB SAMPLEPL=$sPL SAMPLECN=$sCN at least one of these fields has invalid values. "
    echo -e "program=$0 stopped at line=$LINENO.\nReason=$MSG" | mail -s "[Task #${reportticket}]" "$redmine,$email"
    exit 1;
fi

RGparms=$( echo "ID=${sID}:LB=${sLB}:PL=${sPL}:PU=${sPU}:SM=${sSM}:CN=${sCN}" )
rgheader=$( echo -n -e "@RG\t" )$( echo -e "${RGparms}"  | tr ":" "\t" | tr "=" ":" )

######### Alignment: The default settings
if [ ! -d $results ]; then
	mkdir $results
fi

cd $results
mkdir default
cd default

START=$(date +%s)
$bwamemdir/bwa mem -M -t $thr -R  "${rgheader}" $bwa_index $read1 $read2 | $samtoolsdir/samtools view -@ $thr -bS > a.default.0.bam 
END=$(date +%s)
[ -s a.default.0.bam ] && echo "Default alignment successeful!" || exit
alignments=$($samtoolsdir/samtools view -@ $thr -c a.default.0.bam)
if [ "$alignments" -eq 0]; then
	echo			
	echo " Unfortunately, I can NOT process the your request with default parameters" 
	echo
	exit
fi

DIFF=$(( $END - $START ))

echo 
echo "BWA Mem aligned :$alignments: using parameter :default: (*=0=)"  > a.default.0.summary.txt
echo 
echo "Execution time is :$DIFF: seconds" >> a.default.0.summary.txt
echo 
$samtoolsdir/samtools flagstat a.default.0.bam >> a.default.0.summary.txt  # Generating summary statistics


######### Alignment: The combinatorial settings: changing a variable at a time, with the objective of having a sense of how things work


#echo -e "\n\n########################################################################################"
#echo -e "#############                CHECKING PARAMETERS                         ###############"
#echo -e "########################################################################################\n\n"

declare -a parameters=(k r w d c D m W A B O E L U T)
declare -a min=(3 .5 20 20 300 .1 20 0 1 1 1 1 1 1 10)
declare -a step=(3 .5 20 20 300 .1 20 3 2 2 2 2 2 3 10)
declare -a max=(60 4 200 200 10000 1 200 30 20 20 20 20 20 40 80)

cd $results

mkdir ${parameters[@]}

echo The parameters being tested and their ranges are given below:
echo paramters: ${parameters[@]}, 
echo minimum  : ${min[@]} 
echo maximum  : ${max[@]}

pos=0
while [ $pos -lt ${#parameters[@]} ]; do
        par=${parameters[pos]}
	cd $results/$par
        for i in $(seq ${min[pos]} ${step[pos]} ${max[pos]}); do
		START=$(date +%s)	
		$bwamemdir/bwa mem -t $thr -$par "$i" -M -R "${rgheader}" $bwa_index  $read1 $read2 | $samtoolsdir/samtools view -bS -@ $thr > "a.$par.$i.bam" 
		END=$(date +%s)
		DIFF=$(( $END - $START ))
		if [ -s "a.$par.$i.bam" ]; then
			echo "Alignment successeful! with -$par $i" 
		else 
			echo "BWA aligned :0: using default parameter (*=0=)"> "a.$par.$i.summary.txt"
			echo "Execution time is :0: seconds" >> "a.$par.$i.summary.txt"
			continue
		fi
		alignments=$($samtoolsdir/samtools view -@ $thr -c a.$par.$i.bam)
		if [ "$alignments" -eq 0]; then
			echo			
			echo " Unfortunately, I can NOT process the parameter $par = $i with bwa mem" > "a.$par.$i.summary.txt"
			echo
			echo "BWA aligned :0: using default parameter (*=0=)" > "a.$par.$i.summary.txt"
			echo "Execution time is :0: seconds" >> "a.$par.$i.summary.txt"
			continue
		fi

		echo 
		echo "BWA Mem aligned :$alignments: using parameter :$par: =$i" > "a.$par.$i.summary.txt"
		echo 
		echo "Execution time is :$DIFF: seconds" >> "a.$par.$i.summary.txt"
		echo 
		$samtoolsdir/samtools flagstat a.$par.$i.bam >> "a.$par.$i.summary.txt"
	done
	let pos+=1
done


echo -e "\n\n########################################################################################"
echo -e "#############               Summarizing the data :)              ###############"
echo -e "########################################################################################\n\n"

cd $results
echo parameter value Total_reads Total_aligned Mean_MAPQ Time > changing_parameters.txt

ls $results > list
readarray parameters < list
rm list

for par in "${parameters[@]}"; do
        cd $results/$par
        bams="a.$par.*.bam"
        bams=$(echo $bams | tr -d ' ')
        ls $bams|sed 's/.bam//g'|sed 's/a.//g'|cut -d'.' -f2- >range # need to account for cases when range is float
        readarray range < range
        i=0
        while [ $i -lt ${#range[@]} ] ; do
                file="a.$par.${range[i]}.summary.txt"
                summaryfile=$(echo $file | tr -d ' ')   #remove the space introduced from the array variable
                total=$(grep total $summaryfile | cut -d ' ' -f1)
                aligned=$(grep BWA $summaryfile|cut -d':' -f2)
                time=$(grep Execution $summaryfile|cut -d':' -f2)
                file="a.$par.${range[i]}.bam"
                bamfile=$(echo $file | tr -d ' ')
                mapq=$($samtoolsdir/samtools view $bamfile | awk '{sum+=$5} END { print sum/NR}')
                echo $par ${range[i]} $total $aligned $mapq $time >> $results/changing_parameters.txt
                let i+=1
        done
done



echo -e "\n\n########################################################################################"
echo -e "##############                 plotting now for this sample=$SampleName                            ##################"   
echo -e "########################################################################################\n\n"

module load R
Rscript $scriptdir/bwaplotting.R  $results/changing_parameters.txt


echo -e "\n\n########################################################################################"
echo -e "##############                 EXITING NOW                            ##################"	
echo -e "########################################################################################\n\n"
