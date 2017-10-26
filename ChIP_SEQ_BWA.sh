# ----------------------------------------------------------
#!/bin/bash
# request Bourne shell as shell for job
#$ -S /bin/bash
# assume current working directory as paths
#$ -cwd
#$ -N STAR_SE_HUMAN
#
# print start date and time
date

# The job to run if not already
## Quality control with fastqc
ls *.fastq.gz | xargs -P 24 fastqc -t 10

# Trim
for file in `ls *.fastq.gz`
do
newname=`basename $file | sed -e "s/.fastq.gz/.NoAdapt.Trim.fastq.gz/"`
oldname=$(echo ${file} | sed 's/.fastq.gz//')
cat "$oldname"_fastqc/fastqc_data.txt  | grep Over -A 100 | grep 'Illumina\|TruSeq' | grep -P '^[A-Z]' | nl | awk '{print ">" $1 "_adapter\n" $2}' > "$oldname".adapters.fa
java -jar /U3/stefano/src/Trimmomatic-0.35/trimmomatic-0.35.jar SE -threads 14 -phred33 "$oldname".fastq.gz "$newname" ILLUMINACLIP:"$oldname".adapters.fa:2:30:10 SLIDINGWINDOW:4:18 TRAILING:7 LEADING:7 MINLEN:35
echo $oldname
echo $newname
done

# Align with BWA 
for file in `ls *.NoAdapt.Trim.fastq.gz`
do
bamname=`basename $file | sed -e "s/.NoAdapt.Trim.fastq.gz/.bam/"`
bwa mem -t 14 -M /U5/Stefano/Genomes/UCSC/hg38_BWA/hg38.fa $file | samtools view -bS - > "$bamname"
echo $bamname
done

# Mapped + MAPQ > 10 (can go higher if you want)
ls *.bam | parallel --progress --eta -j 15 'samtools view -bq 10 {} > {.}_mapQ10.bam'

# Sort Bam Files 
ls *_mapQ10.bam | parallel --progress --eta -j 15 'samtools sort {} > {.}_sorted.bam'
rename s/_mapQ10_sorted.bam/_sorted.bam/ *

# Remove Dups
for file in `ls *_sorted.bam`
do
newname=`basename $file | sed -e "s/_sorted.bam/.NoDups.bam/"`
matrices=`basename $file | sed -e "s/_sorted.bam/.DupMatrices.txt/"`
java -Xmx4g -jar /U1/zoello/src/picard/picard-tools-1.77/picard-tools-1.77/MarkDuplicates.jar INPUT=$file OUTPUT="$newname" METRICS_FILE="$matrices" VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE TMP_DIR=~/tmp
echo $file
echo $newname
echo $matrices
done

# MACS2 calls
# Make a sample file for you data.
# e.g.
# REP1
# REP2
# REP3

ls *.NoDups.bam | sed "s/_IP.NoDups.bam//;s/_IgG.NoDups.bam//;s/_input.NoDups.bam//;" | awk '!seen[$0]++' > sample_name.txt

## put the unique sample names into an array 
sample_files=($(cut -f 1 sample_name.txt))
# print out all the element in the array
echo "${sample_files[@]}"

## loop over the samples and call peak with macs2  
# For IP
for file in "${sample_files[@]}"
do
	IP_bam="${file}"_IP.NoDups.bam
	Input_bam="${file}"_input.NoDups.bam
	# call regular sharp peaks
	macs2 callpeak -t "$IP_bam" -c "$Input_bam" -g hs -n "${file}"_IP_Model -f BAM -g hs --bw 300 --tsize 76 --extsize 300 -q 0.05 --nolambda
	macs2 callpeak -t "$IP_bam" -c "$Input_bam" -g hs -n "${file}"_IP_NoModel -f BAM -g hs --nomodel --bw 300 --tsize 76 --extsize 300 -q 0.05 --nolambda
done

# For IgG (if you have any):
for file in "${sample_files[@]}"
do
	IP_bam="${file}"_IgG.NoDups.bam
	Input_bam="${file}"_input.NoDups.bam
	# call regular sharp peaks
	macs2 callpeak -t "$IP_bam" -c "$Input_bam" -g hs -n "${file}"_IgG_Model -f BAM -g hs --bw 300 --tsize 76 --extsize 300 -q 0.05 --nolambda
	macs2 callpeak -t "$IP_bam" -c "$Input_bam" -g hs -n "${file}"_IgG_NoModel -f BAM -g hs --nomodel --bw 300 --tsize 76 --extsize 300 -q 0.05 --nolambda
done

# Remove the IgG peaks from the IP samples
cp *_NoModel/*narrowPeak ./
for file in "${sample_files[@]}"
do
	IP_peaks="${file}"_IP_NoModel_peaks.narrowPeak
	IgG_peaks="${file}"_IgG_NoModel_peaks.narrowPeak
 	intersectBed -a "$IP_peaks" -b "$IgG_peaks" -v > ${i}.exclusivePeaks.bed
	annotatePeaks.pl ${i}.exclusivePeaks.bed hg38 -size "given" > ${i}.anno.txt 
done

# Run ChIP seq Quality Metrices (SPP should be installed)
for file in "${sample_files[@]}"
do
	IP_bam="${file}"_IP.NoDups.bam
	Rscript /U3/stefano/SCRIPTS/run_spp.R -c="$Input_bam" -savp="$Input_bam".spp.pdf -out="$Input_bam".run_spp.out  
done                             

## Track files
# Need hg38.chrom.sizes file and .pl script
fetchChromSizes hg38 > hg38.chrom.sizes
for i in *.NoDups.bam
do 
makeTagDirectory $i"TagDir" $i -keepAll
makeUCSCfile $i"TagDir" -o auto -fsize 1e50 -res 10 -norm 1e7 -strand both
cd $i"TagDir"
gunzip -c $i"TagDir".ucsc.bedGraph.gz > $i.ucsc.bedGraph
rm $i"TagDir".ucsc.bedGraph.gz
perl /U3/stefano/SCRIPTS/removeOutOfBoundsReadsFixed.pl $i.ucsc.bedGraph hg38 -chromSizes ../hg38.chrom.sizes > $i.ucsc.fixed.bedGraph
rm $i.ucsc.bedGraph
bedGraphToBigWig $i.ucsc.fixed.bedGraph ../hg38.chrom.sizes $i.bw
rm $i.ucsc.fixed.bedGraph
cd ..
done

# print end date and time again
date
#------------------------------------





















