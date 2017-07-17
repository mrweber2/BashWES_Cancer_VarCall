#PBS -S /bin/bash 			
#PBS -q default				
#PBS -l nodes=1:ppn=3
#PBS -o /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/src-azza/outputs.log
#PBS -e /home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/src-azza/errors.log
#PBS -M aeahmed@illinois.edu
#PBS -m abe

referencedir="/home/groups/hpcbio_shared/azza/H3A_NextGen_assessment_set3/data/genome"

for i in `seq 1 22` ; do

  /home/apps/java/jdk1.8.0_65/bin/java -jar /home/apps/gatk/gatk-3.6/GenomeAnalysisTK.jar \
    -R ${referencedir}/ucsc.hg19.fasta \
    -T SelectVariants \
    -V ${referencedir}/1000G_phase1.indels.hg19.sites.vcf \
    -L chr${i} \
    -o ${referencedir}/IndelsByChr/1000G.chr${i}.vcf
  echo "splitting the 1000G file;       chr: chr${i};         exit with $?"

  /home/apps/java/jdk1.8.0_65/bin/java -jar /home/apps/gatk/gatk-3.6/GenomeAnalysisTK.jar \
    -R ${referencedir}/ucsc.hg19.fasta \
    -T SelectVariants \
    -V ${referencedir}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
    -L chr${i} \
    -o ${referencedir}/IndelsByChr/Mills.chr${i}.vcf

  echo "splitting the Mills file;       chr: chr${i};         exit with $?"

done


for i in chrX chrY chrM ; do

  /home/apps/java/jdk1.8.0_65/bin/java -jar /home/apps/gatk/gatk-3.6/GenomeAnalysisTK.jar \
    -R ${referencedir}/ucsc.hg19.fasta \
    -T SelectVariants \
    -V ${referencedir}/1000G_phase1.indels.hg19.sites.vcf \
    -L ${i} \
    -o ${referencedir}/IndelsByChr/1000G.${i}.vcf
  echo "splitting the 1000G file;       chr: chr${i};         exit with $?"

  /home/apps/java/jdk1.8.0_65/bin/java -jar /home/apps/gatk/gatk-3.6/GenomeAnalysisTK.jar \
    -R ${referencedir}/ucsc.hg19.fasta \
    -T SelectVariants \
    -V ${referencedir}/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
    -L ${i} \
    -o ${referencedir}/IndelsByChr/Mills.${i}.vcf
  echo "splitting the Mills file;       chr: chr${i};         exit with $?"

done

