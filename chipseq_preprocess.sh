#!/bin/bash/

fq=$1
output_Dir=$2
threads=
genome=


base=`basename $fq | sed s/.fastq.gz//g`
echo "Sample name is $base"    
echo "Output directory is $output_Dir"
echo "Reference is ${genome}"


mkdir -p ${output_Dir}/bowtie2/intermediate_bams
mkdir -p ${output_Dir}/trimmedRead

# set up output filenames and locations

## set up file names
align_out=${output_Dir}/bowtie2/${base}_unsorted.sam
align_bam=${output_Dir}/bowtie2/${base}_unsorted.bam
align_sorted=${output_Dir}/bowtie2/${base}_sorted.bam
align_filtered=${output_Dir}/bowtie2/${base}_rmdup.bam
trimmed_read=${output_Dir}/trimmedRead/${base}_trimmed.fq.gz

## set up more variables for 2 additional directoties to help clean up the results folder
bowtie_results=${output_Dir}/bowtie2
intermediate_bams=${output_Dir}/bowtie2/intermediate_bams

# Run Trim_Galore
echo "Run Trim_Galore"
trim_galore -j 4 --fastqc -o ${output_Dir}/trimmedRead $fq

# Run bowtie2
echo "Running bowtie2"
bowtie2 -p $threads -q -x $genome -U $trimmed_read -S $align_out 2> ${output_Dir}/bowtie2/${base}_align.log

# Create BAM from SAM
samtools view -h -S -b -@ $threads -o $align_bam $align_out

# Sort BAM file by genomic coordinates
sambamba sort -q -t $threads -o $align_sorted $align_bam

# Filter out duplicates
echo "Filtering out duplicates"
sambamba view -q -h -t $threads -f bam -F "[XS] == null and not unmapped and not duplicate" $align_sorted > $align_filtered
echo "Filtered by 'null and not unmapped and not duplicate'"    

# Create indices for all the bam files for visualization and QC
samtools index -@ $threads $align_filtered

# Move intermediate files out of the bowtie2 directory
mv $bowtie_results/${base}*sorted* $intermediate_bams
echo "Finished! Ready for peak calling!"
