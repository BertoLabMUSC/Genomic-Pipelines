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

## Needed tools: 
# FASTQC
# Trimmomatic
# STAR
# samtools
# RseqQC
# htseq-count

## Files you need
# hg38.fa = download from USCS
# gencode.v24.annotation.gtf or similar = from gene code. This is the gene annotation (full). In Linux: wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.annotation.gtf.gz; gzip gencode.v24.annotation.gtf.gz > gencode.v24.annotation.gtf 
# genecode_hg38_protein.coding.gtf = filtered gtf file from gencode.v24.annotation.gtf. To filter use grep "protein_coding" gencode.v24.annotation.gtf > genecode_hg38_protein.coding.gtf
# GRCh38_rRNA.bed = from RseQC website
# hg38_UCSC_knownGene.bed = from RseQC website

## Create STAR index 
# Include the GTF file for a better junction approximation
# Overhandg is approximately your read length (in our case Nextseq = 75).
# Download the genome (hg38.fa) from UCSC and filter for only KNOWN chromosomes.
# You can download each chromosome and after using a simple bash script: cat chr*.fa > hg38.fa 
# Download the gene annotation: gencode.v24.annotation.gtf
# Download STAR aligner (for RNA-seq). 
# STAR --runMode genomeGenerate --genomeDir hg38_STAR/ --genomeFastaFiles hg38_STAR/*.fa --runThreadN 13 --sjdbGTFfile hg38/gencode.v24.annotation.gtf --sjdbOverhang 75

## Quality control with fastqc
ls *.fastq.gz | xargs -P 24 fastqc -t 10

## Trim for quality using Trimmomatic. 
# Note: trimmomatic should be in the path: directory were is located/installed the tool.
for file in `ls *.fastq.gz`
do
newname=`basename $file | sed -e "s/.fastq.gz/.NoAdapt.Trim.fastq.gz/"`
oldname=$(echo ${file} | sed 's/.fastq.gz//')
cat "$oldname"_fastqc/fastqc_data.txt  | grep Over -A 100 | grep 'Illumina\|TruSeq' | grep -P '^[A-Z]' | nl | awk '{print ">" $1 "_adapter\n" $2}' > "$oldname".adapters.fa
java -jar /PATH/Trimmomatic-0.35/trimmomatic-0.35.jar SE -threads 14 -phred33 "$oldname".fastq.gz "$newname" ILLUMINACLIP:"$oldname".adapters.fa:2:30:10 SLIDINGWINDOW:4:18 TRAILING:7 LEADING:7 MINLEN:35
echo $oldname
echo $newname
done

# STAR allowing 3 mismatches
mkdir ALIGNMENT/
for file in `ls *.NoAdapt.Trim.fastq.gz`
do
outputname=`basename $file | sed -e "s/.NoAdapt.Trim.fastq.gz/_OUTPUT/"`
STAR --runThreadN 14 --genomeDir /PATH/Genomes/UCSC/hg38_STAR/ --readFilesIn $file --readFilesCommand zcat --sjdbGTFfile /PATH/Genomes/UCSC/hg38_STAR/gencode.v24.annotation.gtf --outFilterType BySJout --outFilterMultimapNmax 10 --alignSJoverhangMin 10 --alignSJDBoverhangMin 1 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outFilterMismatchNmax 3 --twopassMode Basic --outFileNamePrefix ALIGNMENT/$outputname
echo $outputname
done

mv ALIGNMENT/*.out.bam ./

# Get a list of multi-mapped reads
mkdir MULTIMAPPEDREADS/
for file in `ls *_R1_001_OUTPUTAligned.sortedByCoord.out.bam`
do
outputname=`basename $file | sed -e "s/_OUTPUTAligned.sortedByCoord.out.bam/.MultiMappedReads.txt/"`
samtools view $file | grep -v NH:i:1 | perl -pe 's/AS.+(NH:i:\d+)/\$1/' | cut -f1,10,12 | perl -pe 's/NH:i://' | sort -u -k3,3nr > MULTIMAPPEDREADS/"$outputname"
echo $outputname
done

# Fetch only primary alignment (remove unmapped,chimeric etc etc)
find . -name "*_OUTPUTAligned.sortedByCoord.out.bam" | xargs -n 1 -P 12 -iFILES sh -c 'samtools view -F 256 -b FILES > FILES.PRIMARY.bam;'; 
rename s/_OUTPUTAligned.sortedByCoord.out.bam\// *.PRIMARY.bam
mv *.out.bam ALIGNMENT/

# Split Primary alignment in rRNA (in.bam) and non rRNA (ex.bam).
# You need the GRCh38_rRNA.bed from RseqQC data base
for file in `ls *.PRIMARY.bam`
do
samname=`basename $file | sed -e "s/.PRIMARY.bam//"`
split_bam.py -i $file -r /PATH/Genomes/UCSC/hg38_STAR/GRCh38_rRNA.bed -o "$samname"
echo $samname
done

rm *.junk.bam

# Get uniquely mapped
for file in `ls *.ex.bam`
do
outputname=`basename $file | sed -e "s/.ex.bam/.UNIQUE.bam/"`
(samtools view -H $file; samtools view -F 2308 $file | grep -w 'NH:i:1') | samtools view -bS - > "$outputname"
echo $file
echo $outputname
done

## Counts reads in bam
mkdir BAMCOUNT/
find . -name "*.bam" | xargs -n 1 -P 12 -iFILES sh -c 'samtools view -c FILES > BAMCOUNT/FILES.Count.txt;'; 
rename s/.bam\// BAMCOUNT/*.txt

# Genome Coverage 
# You need the hg38_UCSC_knownGene.bed from RseqQC data base
mkdir COVERAGE/
for file in `ls *.UNIQUE.bam`
do
newname=`basename $file | sed -e "s/.bam//"`
geneBody_coverage.py -i $file -r /PATH/Genomes/UCSC/hg38_STAR/hg38_UCSC_knownGene.bed -o COVERAGE/"$newname"
echo $file
echo $newname
done

# Read Distribution
# You need the hg38_UCSC_knownGene.bed from RseqQC data base
mkdir DISTRIBUTION/
find . -name "*.UNIQUE.bam" | xargs -n 1 -P 4 -iFILES sh -c 'read_distribution.py -i FILES -r /PATH/Genomes/UCSC/hg38_STAR/hg38_UCSC_knownGene.bed > DISTRIBUTION/FILES.DistributionLog.txt;';
rename s/UNIQUE.bam\// DISTRIBUTION/*.txt
  
##HTseq-count
# This script count the reads accoring to the GTF that goes along with the genome. 
# Change the -s option based on the RNA seq data you have
# Change the -t option based on the RNA seq data (polyA = exon, Total RNA = gene)
mkdir HTSEQ/
parallel -j 14 'samtools view {} | htseq-count -m intersection-strict -t exon -i gene_name -s reverse - /PATH/Genomes/UCSC/hg38_STAR/genecode_hg38_protein.coding.gtf > {.}.txt' ::: *.UNIQUE.bam

# Remove list lines from HTseq count
for file in `ls *.txt`
do
newname=`basename $file | sed -e "s/.txt/.LastLinesRem.txt/"`
head -n -5 $file > "$newname"
echo $file
echo $newname
done

mv *.txt HTSEQ/

# print end date and time again
date
#------------------------------------

