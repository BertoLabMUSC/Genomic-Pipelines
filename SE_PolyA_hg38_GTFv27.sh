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
# gencode.v27.annotation.gtf or similar = from gene code. This is the gene annotation (full). In Linux: wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz; gzip gencode.v27.annotation.gtf.gz > gencode.v27.annotation.gtf 
# gencode.v27.protein_coding.gtf = filtered gtf file from gencode.v24.annotation.gtf. To filter use grep "protein_coding" gencode.v27.annotation.gtf > gencode.v27.protein_coding.gtf
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

# Paths to the genome directories
DBDIR=/U5/Stefano/Genomes/UCSC/hg38_STAR_GencodeV27/
RSEMDIR=/U5/Stefano/Genomes/UCSC/hg38_RSEM_GencodeV27/hg38

## Quality control with fastqc
ls *.fastq.gz | xargs -P 24 fastqc -t 10

## Trim for quality using Trimmomatic. 
# Note: trimmomatic should be in the path: directory were is located/installed the tool.
for file in `ls *.fastq.gz`
do
newname=`basename $file | sed -e "s/.fastq.gz/.NoAdapt.Trim.fastq.gz/"`
oldname=$(echo ${file} | sed 's/.fastq.gz//')
cat "$oldname"_fastqc/fastqc_data.txt  | grep Over -A 100 | grep 'Illumina\|TruSeq' | grep -P '^[A-Z]' | nl | awk '{print ">" $1 "_adapter\n" $2}' > "$oldname".adapters.fa
java -jar /U3/stefano/src/Trimmomatic-0.36/Trimmomatic-0.36.jar SE -threads 14 -phred33 "$oldname".fastq.gz "$newname" ILLUMINACLIP:"$oldname".adapters.fa:2:30:10 SLIDINGWINDOW:4:18 TRAILING:7 LEADING:7 MINLEN:35
echo $oldname
echo $newname
done

# STAR allowing 3 mismatches
mkdir ALIGNMENT/
for file in `ls *.NoAdapt.Trim.fastq.gz`
do
outputname=`basename $file | sed -e "s/.NoAdapt.Trim.fastq.gz//"`
STAR --runThreadN 14 \
--genomeDir /U5/Stefano/Genomes/UCSC/hg38_STAR_GencodeV27/ \
--readFilesIn $file \
--readFilesCommand zcat \
--sjdbGTFfile /U5/Stefano/Genomes/UCSC/hg38_STAR_GencodeV27/gencode.v27.annotation.gtf \
--outFilterType BySJout  \
--outFilterMismatchNoverReadLmax 0.04 \
--outFilterMultimapNmax 10 \
--alignSJoverhangMin 10 \
--alignSJDBoverhangMin 1 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outFilterMismatchNmax 3 \
--twopassMode Basic \
--outFileNamePrefix ALIGNMENT/$outputname \
--chimSegmentMin 15 \
--chimScoreMin 15 \
--chimScoreSeparation 10 \
--chimJunctionOverhangMin 15 \
--quantMode TranscriptomeSAM
echo $outputname
done

mv ALIGNMENT/*.out.bam ./

# Fetch only primary alignment (remove unmapped,chimeric etc etc)
find . -name "*Aligned.sortedByCoord.out.bam" | xargs -n 1 -P 12 -iFILES sh -c 'samtools view -F 256 -b FILES > FILES.PRIMARY.bam;'; 
rename s/Aligned.sortedByCoord.out.bam\// *.PRIMARY.bam
mv *.out.bam ALIGNMENT/

# Split Primary alignment in rRNA (in.bam) and non rRNA (ex.bam).
# You need the GRCh38_rRNA.bed from RseqQC data base
for file in `ls *.PRIMARY.bam`
do
samname=`basename $file | sed -e "s/.PRIMARY.bam//"`
split_bam.py -i $file -r /U5/Stefano/Genomes/UCSC/hg38_STAR_GencodeV27/GRCh38_rRNA.bed -o "$samname"
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
geneBody_coverage.py -i $file -r /U5/Stefano/Genomes/UCSC/hg38_STAR_GencodeV27/hg38_UCSC_knownGene.bed -o COVERAGE/"$newname"
echo $file
echo $newname
done

# Read Distribution
# You need the hg38_UCSC_knownGene.bed from RseqQC data base
mkdir DISTRIBUTION/
find . -name "*.UNIQUE.bam" | xargs -n 1 -P 4 -iFILES sh -c 'read_distribution.py -i FILES -r /U5/Stefano/Genomes/UCSC/hg38_STAR_GencodeV27/hg38_UCSC_knownGene.bed > DISTRIBUTION/FILES.DistributionLog.txt;';
rename s/UNIQUE.bam\// DISTRIBUTION/*.txt
  
##HTseq-count
# This script count the reads accoring to the GTF that goes along with the genome. 
# Change the -s option based on the RNA seq data you have
# Change the -t option based on the RNA seq data (polyA = exon, Total RNA = gene)
mkdir HTSEQ/
parallel -j 14 'samtools view {} | htseq-count -m intersection-strict -t exon -i gene_name -s reverse - /U5/Stefano/Genomes/UCSC/hg38_STAR_GencodeV27/gencode.v27.protein_coding.gtf > {.}.txt' ::: *.UNIQUE.bam

# Remove list lines from HTseq count
for file in `ls *.txt`
do
newname=`basename $file | sed -e "s/.txt/.LastLinesRem.txt/"`
head -n -5 $file > "$newname"
echo $file
echo $newname
done
mv *.txt HTSEQ/

# RSEM quantification
mv ALIGNMENT/*.toTranscriptome.out.bam ./
mkdir RSEM/
for file in `ls *Aligned.toTranscriptome.out.bam`
do
newname=`basename $file | sed -e "s/Aligned.toTranscriptome.out.bam/.rsem/"`
rsem-calculate-expression --bam --estimate-rspd  --calc-ci --no-bam-output --seed 12345 -p 14 --ci-memory 30000 --forward-prob 0 $file /U5/Stefano/Genomes/UCSC/hg38_RSEM_GencodeV27/hg38 RSEM/"$newname" 
echo $file
echo $newname
done
mv *.toTranscriptome.out.bam ALIGNMENT/

## Track files both strand
# Need hg38.chrom.sizes file and .pl script
fetchChromSizes hg38 > hg38.chrom.sizes
for i in *.UNIQUE.bam
do 
makeTagDirectory $i"TagDir" $i -keepAll
makeUCSCfile $i"TagDir" -o $i"TagDir"/$i'.POS.bedGraph' -fsize 1e50 -res 10 -norm 1e7 -strand -
makeUCSCfile $i"TagDir" -o $i"TagDir"/$i'.NEG.bedGraph' -fsize 1e50 -res 10 -norm 1e7 -strand +
cd $i"TagDir"
gunzip -c $i'.NEG.bedGraph'.gz | awk '{printf "%s\t%s\t%s\t-%s\n", $1, $2, $3, $4}' > $i'.NEG.bedGraph'
gunzip -c $i'.POS.bedGraph'.gz | awk '{printf "%s\t%s\t%s\t+%s\n", $1, $2, $3, $4}' > $i'.POS.bedGraph' 
perl /U3/stefano/SCRIPTS/removeOutOfBoundsReadsFixed.pl $i'.POS.bedGraph' hg38 -chromSizes ../hg38.chrom.sizes > $i'.POS.fix.bedGraph'
perl /U3/stefano/SCRIPTS/removeOutOfBoundsReadsFixed.pl $i'.NEG.bedGraph' hg38 -chromSizes ../hg38.chrom.sizes > $i'.NEG.fix.bedGraph'
bedGraphToBigWig $i'.POS.fix.bedGraph' ../hg38.chrom.sizes $i.for.bw
bedGraphToBigWig $i'.NEG.fix.bedGraph' ../hg38.chrom.sizes $i.rev.bw
rm *.bedGraph
rm *.bedGraph.gz
cd ..
done

rm *.ex.bam
rm *.in.bam
mkdir TRACKS/
mv *TagDir* TRACKS/
mkdir FASTQ/
mv *.fastq.gz FASTQ/
mkdir FASTQC/
mv *_fastqc* FASTQC/
mkdir BAMS/
mv *.bam BAMS/

# print end date and time again
date
#--------------------