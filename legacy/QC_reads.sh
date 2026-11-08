#!/bin/bash

#NOTES: Feed positional arguments as noted below.

unset genome_R1
unset genome_R2
unset numthreads
unset trimmomaticFolder
unset trimmomaticJar
unset outdirectory
unset outname
unset flashProg
unset pigzProg
unset qcset

numthreads=$1 #Number of threads to use
outname=$2 #Simple name for appending to output files
outdirectory=$3 #Location where files should be stored (will be made at start)
genome_R1=$4 #Input raw reads from first read pair (R1) (absolute path)
genome_R2=$5 #Input raw reads from second read pair (R2) (absolute path)
trimmomaticJar=$6 #trimmomatic in path
trimmomaticFolder=$7 #Path to the trimmomatic folder, which contains adapters/TruSeq3-PE-2.fa
flashProg=$8 #Full path to flash binary
pigzProg=$9 #Full path to pigz binary
qcset=${10} #QC level

mkdir -p $outdirectory
touch ${outdirectory}/${outname}_qc_run.log
myInvocation="$(printf %q "$BASH_SOURCE")$((($#)) && printf ' %q' "$@")"

exec &> >(tee -a ${outdirectory}/${outname}_qc_run.log)
echo
echo "Start of run:"
date
echo
echo "Invoked script options:"
echo "$myInvocation"
echo

##########################################################################################
##
##    QUALITY CONTROL READS (Trimmomatic & Flash)
##
##########################################################################################

#Remove adaptors
$trimmomaticJar PE -threads $numthreads -baseout $outdirectory/$outname $genome_R1 $genome_R2 ILLUMINACLIP:${trimmomaticFolder}/adapters/TruSeq3-PE-2.fa:2:30:10 >> ${outdirectory}/${outname}_qc_run.log 2>&1

#Clean unpaired reads
$trimmomaticJar SE -threads $numthreads $outdirectory/${outname}_1U $outdirectory/${outname}_1U.qual_trimmed SLIDINGWINDOW:4:${qcset} MINLEN:50 >> $outdirectory/${outname}_qc_run.log 2>&1

$trimmomaticJar SE -threads $numthreads $outdirectory/${outname}_2U $outdirectory/${outname}_2U.qual_trimmed SLIDINGWINDOW:4:${qcset} MINLEN:50 >> $outdirectory/${outname}_qc_run.log 2>&1

#Merge paired reads
$flashProg --threads $numthreads --output-prefix $outname --max-overlap 150 $outdirectory/${outname}_1P $outdirectory/${outname}_2P --output-directory $outdirectory >> $outdirectory/${outname}_qc_run.log 2>&1

#mv flashed_reads/flashed.extendedFrags.fastq reads.adapter_trimmed.flash_combined
#mv flashed_reads/flashed.notCombined_1.fastq reads.adapter_trimmed.flash_notcombined_1P
#mv flashed_reads/flashed.notCombined_2.fastq reads.adapter_trimmed.flash_notcombined_2P

#Clean merged reads
$trimmomaticJar SE -threads $numthreads $outdirectory/$outname.extendedFrags.fastq $outdirectory/$outname.extendedFrags.final.fastq SLIDINGWINDOW:4:${qcset} MINLEN:50 >> ${outname}_qc_run.log 2>&1

#Clean unmerged, paired reads
$trimmomaticJar PE -threads $numthreads $outdirectory/$outname.notCombined_1.fastq $outdirectory/$outname.notCombined_2.fastq -baseout  $outdirectory/$outname.notCombined.final SLIDINGWINDOW:4:${qcset} MINLEN:50 >> ${outname}_qc_run.log 2>&1

##move files into ${outname}_1.fq, ${outname}_2.fq, and ${outname}_unpaird.fq files for downstream use
#mv reads.adapter_trimmed.flash_notcombined.qual_trimmed_1P ${outname}_reads.1.fq
#mv reads.adapter_trimmed.flash_notcombined.qual_trimmed_2P ${outname}_reads.2.fq

cat $outdirectory/$outname.extendedFrags.final.fastq $outdirectory/$outname.notCombined.final_1U $outdirectory/$outname.notCombined.final_2U $outdirectory/${outname}_1U.qual_trimmed $outdirectory/${outname}_2U.qual_trimmed > $outdirectory/${outname}_reads.U.fq

#GZIP files (pigz preferred)
$pigzProg --best --processes $numthreads $outdirectory/$outname.notCombined.final_1P $outdirectory/$outname.notCombined.final_2P $outdirectory/${outname}_reads.U.fq

#mkdir temp_delete
#mv reads.adapter_trimmed* temp_delete/
#rm -rf temp_delete/
mkdir $outdirectory/cleaned
mv $outdirectory/*gz $outdirectory/cleaned/
rm $outdirectory/${outname}*
mv $outdirectory/cleaned/*gz $outdirectory/
#STATS File?

sleep 1
