#!/bin/bash
# Created : 2025. 03.
# Script by Chanyeon Park.
# Last Update : 2025. 03.


while getopts ":s:t:r:y:p:o:" opt; do
  case $opt in
    s) INPUT="$OPTARG" ;;
    t) THREAD="$OPTARG" ;;
    r) GENOME_VERSION="$OPTARG" ;;
    y) SEQ_TYPE="$OPTARG" ;;
    p) RUNNING_STEP="$OPTARG" ;;
    o) EXP_DIR="$OPTARG" ;;
    \?) echo "X Invalid option: -$OPTARG" >&2; exit 1 ;;
    :) echo "X Option -$OPTARG requires an argument." >&2; exit 1 ;;
  esac
done

# validate essential parameters
if [ -z "$INPUT" ] || [ -z "$GENOME_VERSION" ] || [ -z "$EXP_DIR" ]; then
    echo "Chanyeon's RNA-seq pipeline"
    echo "Usage: $(basename $0) -s [sample name] -r [genome version] -o [output directory] [-t threads] [-y sequencing type] [-p step]"
    echo "Example: $(basename $0) -s FMD_1 -t 8 -r mm10 -o RNA_res"
    echo "    -s    If you have the FMD_1.fastq.gz file, sample name is FMD_1 only. Not FMD_1.fastq.gz."
    echo "    -r    mm10, hg19, hg38 references are currently supported."
    echo "    -y    Sequencing type (default: paired): single or paired"
    echo "    -p    Example of step: all (default), trimming, mapping, homer, visualization, postQC, geneCount"
    echo "    If you want to run multiple steps, concatenate step names without spaces (e.g., mapping,homer,visualization,postQC)"
    exit 1
fi

if [ -z "$SEQ_TYPE" ]; then
    SEQ_TYPE="paired"
fi
if [ -z "$RUNNING_STEP" ]; then
    RUNNING_STEP="all"
fi
if [ -z "$THREAD" ]; then
    THREAD=1
fi

SAMPLE_ID=`basename "$INPUT"`
SAMPLE_PATH=`dirname "$INPUT"`



echo "============================================================================================"
echo "         RNA-seq Analysis Pipeline           "
echo "============================================================================================"
echo "   Parameters detected:"
echo "   Sample ID     : $SAMPLE_ID"
echo "   Sample path   : $SAMPLE_PATH"
echo "   Threads       : $THREAD"
echo "   Reference     : $GENOME_VERSION"
echo "   Output dir    : $EXP_DIR"
echo "   Sequencing    : $SEQ_TYPE (default: paired if not specified)"
echo "   Step          : $RUNNING_STEP"
echo "============================================================================================"

#############################################################################
########################### Configuration (Tool)  ###########################
#############################################################################

#############################################################################
######################### Configuration (Reference) #########################
#############################################################################

if [ "$GENOME_VERSION" = "mm10" ]
then
	GTF_FILE=/home/tinyimp/Reference/GRCm38/gencode.vM23.annotation.gtf
	REF_DIR=/home/tinyimp/Reference/GRCm38/STAR
	BED_FILE=/home/tinyimp/Reference/GRCm38/mm10_GENCODE_VM25_basic.bed
	GENOME_GRCm38="ID"
elif [ "$GENOME_VERSION" = "hg19" ]
then
	GTF_FILE=/home/tinyimp/Reference/GRCh37/gencode.v19.annotation.gtf
	REF_DIR=/home/tinyimp/Reference/GRCh37/RNA_reference/
 	BED_FILE=/home/tinyimp/Reference/GRCh37/hg19_GENCODE_V42_Basic.bed
	GENOME_ID="GRCh37"
elif [ "$GENOME_VERSION" = "hg38" ]
then
	GTF_FILE=/home/tinyimp/Reference/GRCh38/gencode.v42.annotation.gtf
	REF_DIR=/home/tinyimp/Reference/GRCh38/STAR/
	BED_FILE=/home/tinyimp/Reference/GRCh38/hg38_GENCODE_V42_Basic.bed
	GENOME_ID="GRCh38"
else
	echo 'Reference Genome Not Detected for Processing. Check the script you used.'
	exit
fi

if [[ -d "$REF_DIR" ]]
then
	echo 'STAR and RSEM reference directory exist. '
else
	echo 'STAR and RSEM reference directory not found. '
	echo 'Check the reference directory. '
	echo 'If proper reference directory not exist, build new reference by embedded script inside the RNA reference directory. '
	exit
fi

#############################################################################
################################ Preparation ################################
#############################################################################

mkdir -p $EXP_DIR/1.RAWDATA
mkdir -p $EXP_DIR/2.QC_report
mkdir -p $EXP_DIR/3.TRIMMED
mkdir -p $EXP_DIR/4.BAM
mkdir -p $EXP_DIR/5.Homer_tag
mkdir -p $EXP_DIR/6.UCSC_file
mkdir -p $EXP_DIR/7.Post_QC
mkdir -p $EXP_DIR/8.RSEM_geneCount

FILE_NUMBER=`ls $SAMPLE_PATH | grep fastq.gz | wc -l`

if [ $FILE_NUMBER -ne 0 ]
then
	mv -n "$SAMPLE_PATH"/"$SAMPLE_ID"*.f*q.gz $EXP_DIR/1.RAWDATA/ 2> /dev/null
fi

#############################################################################
######################### Paired-end or Single-end ##########################
#############################################################################

FILE_NUMBER=`ls $EXP_DIR/1.RAWDATA | grep "$SAMPLE_ID" | grep gz | wc -l`

if [ $FILE_NUMBER = 1 ]
then
	MODE=single
elif [ $FILE_NUMBER = 2 ]
then
	MODE=paired
else
	echo "Check Your Read File List"
fi

if [ "$SEQ_TYPE" = "$MODE" ]
then
	echo "$MODE end detected"
else
	echo "Sequencing type is not matched to given file number"
	echo "Sequencing type: $SEQ_TYPE / Raw data type: $MODE"
	echo "Check Your Read File List"
	exit
fi

#############################################################################
######################## Running step determination #########################
#############################################################################

if [ "$RUNNING_STEP" = "" -o "$RUNNING_STEP" = "all" ]
then
	RUNNING_STEP_TRIM=1
	RUNNING_STEP_STAR=1
	#RUNNING_STEP_STRAND_SPLIT=1
	RUNNING_STEP_HOMER=1
	RUNNING_STEP_VISUALIZE=1
	RUNNING_STEP_POSTQC=1
	RUNNING_STEP_RSEM=1
else
	RUNNING_STEP_TRIM=0
	RUNNING_STEP_STAR=0
	#RUNNING_STEP_STRAND_SPLIT=0
	RUNNING_STEP_HOMER=0
	RUNNING_STEP_VISUALIZE=0
	RUNNING_STEP_POSTQC=0
	RUNNING_STEP_RSEM=0

	RUNNING_STEP_TRIM=`echo $RUNNING_STEP | grep trimming | wc -l`
	RUNNING_STEP_STAR=`echo $RUNNING_STEP | grep mapping | wc -l`
	#RUNNING_STEP_STRAND_SPLIT=`echo $RUNNING_STEP | grep strandSplitting | wc -l`
	RUNNING_STEP_HOMER=`echo $RUNNING_STEP | grep homer | wc -l`
	RUNNING_STEP_VISUALIZE=`echo $RUNNING_STEP | grep visualization | wc -l`
	RUNNING_STEP_POSTQC=`echo $RUNNING_STEP | grep postQC | wc -l`
	RUNNING_STEP_RSEM=`echo $RUNNING_STEP | grep geneCount | wc -l`
fi

#############################################################################
#############################################################################
#############################################################################

if [ $RUNNING_STEP_TRIM = 1 ]
then
	echo "Step1: fastqc and trimming"


#############################################################################
############################## Pre Processing ###############################
#############################################################################

#fastQC for raw reads

fastqc \
        -t $THREAD \
        $EXP_DIR/1.RAWDATA/"$SAMPLE_ID"*.f*q.gz \
        -o \
        $EXP_DIR/2.QC_report

#Trimming

if [ "$MODE" = "single" ]
then
        trim_galore -j 4 \
        --fastqc \
        $EXP_DIR/1.RAWDATA/"$SAMPLE_ID".f*q.gz \
        -o \
        $EXP_DIR/3.TRIMMED/
elif [ "$MODE" = "paired" ]
then
        trim_galore -j 4 \
        --fastqc \
                --paired \
                $EXP_DIR/1.RAWDATA/"$SAMPLE_ID"_1.f*q.gz \
                $EXP_DIR/1.RAWDATA/"$SAMPLE_ID"_2.f*q.gz \
                -o \
                $EXP_DIR/3.TRIMMED/
fi


else
        echo 'Step1: fastqc and trimming (skip)'
fi

#############################################################################
#############################################################################
#############################################################################

if [ $RUNNING_STEP_STAR = 1 ]
then
	echo 'Step2: mapping'

#############################################################################
############################ Mapping to Genome ##############################
#############################################################################

##Check whether You did Generate Genome indices, STAR need to make Genome indices.  

echo "Started STAR for "$GENOME_VERSION" reference"

if [ "$MODE" = "single" ]
then
	echo 'single end mapping'
	STAR \
    --runMode alignReads \
    --alignMatesGapMax 1000000 \
		--runThreadN $THREAD \
		--genomeDir $REF_DIR/ \
		--sjdbGTFfile $GTF_FILE \
		--readFilesCommand zcat \
		--readFilesIn $EXP_DIR/3.TRIMMED/"$SAMPLE_ID"*.f*q.gz \
		--outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM \
		--outFileNamePrefix $EXP_DIR/4.BAM/"$SAMPLE_ID"."$GENOME_VERSION"_ \


elif [ "$MODE" = "paired" ]
then
	echo 'paired end mapping'
	STAR \
    --runMode alignReads \
    --alignMatesGapMax 1000000 \
		--runThreadN $THREAD \
    --genomeDir $REF_DIR/ \
    --sjdbGTFfile $GTF_FILE \
    --readFilesCommand zcat \
    --readFilesIn $EXP_DIR/3.TRIMMED/"$SAMPLE_ID"*_1.f*q.gz $EXP_DIR/3.TRIMMED/"$SAMPLE_ID"*_2.f*q.gz \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode TranscriptomeSAM \
    --outFileNamePrefix $EXP_DIR/4.BAM/"$SAMPLE_ID"."$GENOME_VERSION"_ \

fi

# Indexing Primary Mapped Reads
samtools index \
	-@ $THREAD \
	$EXP_DIR/4.BAM/"$SAMPLE_ID"."$GENOME_VERSION"_Aligned.sortedByCoord.out.bam

infer_experiment.py -r "$BED_FILE" -i $EXP_DIR/4.BAM/"$SAMPLE_ID"."$GENOME_VERSION"_Aligned.sortedByCoord.out.bam > $EXP_DIR/4.BAM/"$SAMPLE_ID"."$GENOME_VERSION"_strandedness_infer.txt

fraction_1=$(grep 'Fraction of reads explained by "1++,1--,2+-,2-+"' "$EXP_DIR"/4.BAM/"$SAMPLE_ID"."$GENOME_VERSION"_strandedness_infer.txt | awk '{print $NF}')
fraction_2=$(grep 'Fraction of reads explained by "1+-,1-+,2++,2--"' "$EXP_DIR"/4.BAM/"$SAMPLE_ID"."$GENOME_VERSION"_strandedness_infer.txt | awk '{print $NF}')

if (( $(echo "$fraction_1 > $fraction_2" | bc -l) )); then
    strandedness="forward"
else
    strandedness="reverse"
fi

echo "Determined strandedness: $strandedness"


else
	echo 'Step2: mapping (skip)'
fi


#############################################################################
#############################################################################
#############################################################################

if [ $RUNNING_STEP_HOMER = 1 ]
then
	echo 'Step3: HOMER tag directory generation'

#############################################################################
####################### Making Homer Tag Directory ##########################
#############################################################################

makeTagDirectory \
	$EXP_DIR/5.Homer_tag/"$SAMPLE_ID"."$GENOME_VERSION" \
	$EXP_DIR/4.BAM/"$SAMPLE_ID"."$GENOME_VERSION"_Aligned.sortedByCoord.out.bam

else
	echo 'Step3: HOMER tag directory generation (skip)'
fi

#############################################################################
#############################################################################
#############################################################################

if [ $RUNNING_STEP_VISUALIZE = 1 ]
then
	echo 'Step4: IGV signal track generation'

#############################################################################
#########################  Making IGV signal Track ##########################
#############################################################################

# Making BigWig file using "deeptools-bamcoverage" , Normalized using CPM

bamCoverage \
	-b $EXP_DIR/4.BAM/"$SAMPLE_ID"."$GENOME_VERSION"_Aligned.sortedByCoord.out.bam \
	-o $EXP_DIR/6.UCSC_file/"$SAMPLE_ID"."$GENOME_VERSION"_Aligned.sortedByCoord.out.CPM.bw \
	--normalizeUsing CPM \
	--exactScaling \
	-p $THREAD

else
	echo 'Step4: IGV signal track generation (skip)'
fi

#############################################################################
#############################################################################
#############################################################################

if [ $RUNNING_STEP_POSTQC = 1 ]
then
	echo 'Step5: Post QC'

#############################################################################
###################### Making Statistical. summary ##########################
#############################################################################

# Mappability statistics

if [ "$MODE" = "single" ]
then
	total_reads_fastq_R1=`zcat $EXP_DIR/1.RAWDATA/"$SAMPLE_ID".fastq.gz | wc -l`
	total_reads_fastq=`expr $total_reads_fastq_R1 / 4`
	echo -e "sample ID\ttotal reads (FASTQ)\ttotal reads (STAR)\tuniquely mapped reads\tmulti-mapped reads\tunique mapping rate\tmulti-mapped mapping rate" > $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION".mappability.stat

elif [ "$MODE" = "paired" ]
then
	total_reads_fastq_R1=`zcat $EXP_DIR/1.RAWDATA/"$SAMPLE_ID"_1.fastq.gz | wc -l`
	total_reads_fastq_R2=`zcat $EXP_DIR/1.RAWDATA/"$SAMPLE_ID"_2.fastq.gz | wc -l`
	

	if [ $total_reads_fastq_R1 -ne $total_reads_fastq_R2 ]
	then
		echo ' FASTQ R1 and R2 are not matched to each other '
		echo ' All of previous results might be wrong. Check your FASTQ file first. '
		exit
	else
		echo -e "sample ID\ttotal reads (FASTQ)\ttotal read pairs (STAR)\tuniquely mapped read pairs\tmulti-mapped read pairs\tunique mapping rate\tmulti-mapped mapping rate" > $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION".mappability.stat
		total_reads_fastq=`expr $total_reads_fastq_R1 / 2`
	fi
fi

total_reads=`cat $EXP_DIR/4.BAM/"$SAMPLE_ID"."$GENOME_VERSION"_Log.final.out | grep Number\ of\ input\ reads | rev | cut -f1 | rev`

mapped_reads_uniq=`cat $EXP_DIR/4.BAM/"$SAMPLE_ID"."$GENOME_VERSION"_Log.final.out | grep Uniquely\ mapped\ reads\ number | rev | cut -f1 | rev`

mapped_reads_multi=`cat $EXP_DIR/4.BAM/"$SAMPLE_ID"."$GENOME_VERSION"_Log.final.out | grep Number\ of\ reads\ mapped\ to\ multiple\ loci | rev | cut -f1 | rev`

mapping_rate_uniq=`cat $EXP_DIR/4.BAM/"$SAMPLE_ID"."$GENOME_VERSION"_Log.final.out | grep Uniquely\ mapped\ reads\ \% | rev | cut -f1 | rev`

mapping_rate_multi=`cat $EXP_DIR/4.BAM/"$SAMPLE_ID"."$GENOME_VERSION"_Log.final.out | grep \%\ of\ reads\ mapped\ to\ multiple\ loci | rev | cut -f1 | rev`

echo -e "$SAMPLE_ID\t$total_reads_fastq\t$total_reads\t$mapped_reads_uniq\t$mapped_reads_multi\t$mapping_rate_uniq\t$mapping_rate_multi" >> $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION".mappability.stat

#############################################################################
############################### Alignment QC ################################
#############################################################################

# For qualimap QC criteria, see https://hbctraining.github.io/Intro-to-rnaseq-fasrc-salmon-flipped/lessons/10_QC_Qualimap.html

qualimap rnaseq \
	--algorithm uniquely-mapped-reads \
	-bam $EXP_DIR/4.BAM/"$SAMPLE_ID"."$GENOME_VERSION"_Aligned.sortedByCoord.out.bam \
	-gtf $GTF_FILE \
	-outfile "$SAMPLE_ID"."$GENOME_VERSION".qualimap \
	-outformat PDF:HTML \
	--sequencing-protocol strand-specific-reverse \
	--paired \
	--java-mem-size="$THREAD"G \
	-outdir $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION"/ 

mv $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION"/rnaseq_qc_results.txt $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION"/"$SAMPLE_ID"."$GENOME_VERSION".qualimap.txt
mv $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION"/qualimapReport.html $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION"/"$SAMPLE_ID"."$GENOME_VERSION".qualimap.html

if [ "$MODE" = "single" ]
then
	qualimap_aligned_reads=`cat $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION"/"$SAMPLE_ID"."$GENOME_VERSION".qualimap.txt | grep total\ alignments | rev | cut -f1 -d " " | rev`
	echo -e "aligned reads (qualimap)\texonic reads\tintronic reads\tintergenic reads\texonic rate\tintronic rate\tintergenic rate" > $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION".qualimap.genomicOrigin.stat

elif [ "$MODE" = "paired" ]
then
	qualimap_aligned_reads=`cat $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION"/"$SAMPLE_ID"."$GENOME_VERSION".qualimap.txt | grep read\ pairs\ aligned | rev | cut -f1 -d " " | rev`
	echo -e "aligned read pairs (qualimap)\texonic read pairs\tintronic read pairs\tintergenic read pairs\texonic rate\tintronic rate\tintergenic rate" > $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION".qualimap.genomicOrigin.stat
fi

qualimap_exonic_reads=`cat $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION"/"$SAMPLE_ID"."$GENOME_VERSION".qualimap.txt | grep exonic | rev | cut -f2 -d " " | rev | sed 's/\([0-9]\),\([0-9]\)/\1\2/g'`
qualimap_exonic_overlap_reads=`cat $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION"/"$SAMPLE_ID"."$GENOME_VERSION".qualimap.txt | grep overlapping\ exon | rev | cut -f2 -d " " | rev | sed 's/\([0-9]\),\([0-9]\)/\1\2/g'`
qualimap_exonic_merged_reads=`expr $qualimap_exonic_reads + $qualimap_exonic_overlap_reads`

qualimap_intronic_reads=`cat $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION"/"$SAMPLE_ID"."$GENOME_VERSION".qualimap.txt | grep intronic | rev | cut -f2 -d " " | rev | sed 's/\([0-9]\),\([0-9]\)/\1\2/g'`
qualimap_intergenic_reads=`cat $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION"/"$SAMPLE_ID"."$GENOME_VERSION".qualimap.txt | grep intergenic | rev | cut -f2 -d " " | rev | sed 's/\([0-9]\),\([0-9]\)/\1\2/g'`

qualimap_exonic_rate=`cat $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION"/"$SAMPLE_ID"."$GENOME_VERSION".qualimap.txt | grep exonic | cut -f2 -d "(" | cut -f1 -d "%"`
qualimap_exonic_overlap_rate=`cat $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION"/"$SAMPLE_ID"."$GENOME_VERSION".qualimap.txt | grep overlapping\ exon |  cut -f2 -d "(" | cut -f1 -d "%"`
qualimap_exonic_merged_rate=`echo "scale=2; $qualimap_exonic_rate + $qualimap_exonic_overlap_rate" | bc`

qualimap_intronic_rate=`cat $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION"/"$SAMPLE_ID"."$GENOME_VERSION".qualimap.txt | grep intronic | cut -f2 -d "(" | cut -f1 -d ")"`
qualimap_intergenic_rate=`cat $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION"/"$SAMPLE_ID"."$GENOME_VERSION".qualimap.txt | grep intergenic | cut -f2 -d "(" | cut -f1 -d ")"`

echo -e "$qualimap_aligned_reads\t$qualimap_exonic_merged_reads\t$qualimap_intronic_reads\t$qualimap_intergenic_reads\t$qualimap_exonic_merged_rate%\t$qualimap_intronic_rate\t$qualimap_intergenic_rate" >> $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION".qualimap.genomicOrigin.stat

# Merge QC files
paste \
	$EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION".mappability.stat \
	$EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION".qualimap.genomicOrigin.stat \
> $EXP_DIR/7.Post_QC/"$SAMPLE_ID"."$GENOME_VERSION".QCvalue.txt

else
	echo 'Step5: Post QC (skip)'
fi

#############################################################################
#############################################################################
#############################################################################

if [ $RUNNING_STEP_RSEM = 1 ]
then
	echo 'Step6: RSEM gene count '

#############################################################################
############################## RSEM Gene Count ##############################
#############################################################################
fraction_1=$(grep 'Fraction of reads explained by "1++,1--,2+-,2-+"' "$EXP_DIR"/4.BAM/"$SAMPLE_ID"."$GENOME_VERSION"_strandedness_infer.txt | awk '{print $NF}')
fraction_2=$(grep 'Fraction of reads explained by "1+-,1-+,2++,2--"' "$EXP_DIR"/4.BAM/"$SAMPLE_ID"."$GENOME_VERSION"_strandedness_infer.txt | awk '{print $NF}')

if (( $(echo "$fraction_1 > $fraction_2" | bc -l) )); then
    strandedness="forward"
else
    strandedness="reverse"
fi



if [ "$MODE" = "single" ]
then
	rsem-calculate-expression \
		--time \
		--no-bam-output \
		-p "$THREAD" \
	  --bam \
    --estimate-rspd \
		$EXP_DIR/4.BAM/"$SAMPLE_ID"."$GENOME_VERSION"_Aligned.toTranscriptome.out.bam \
		$REF_DIR/"$GENOME_VERSION" \
		$EXP_DIR/8.RSEM_geneCount/"$SAMPLE_ID"."$GENOME_VERSION" \
		2> $EXP_DIR/8.RSEM_geneCount/"$SAMPLE_ID"."$GENOME_VERSION".rsem.log

elif [ "$MODE" = "paired" ]
then
	rsem-calculate-expression \
		--strandedness "$strandedness" \
		--time \
		--no-bam-output \
		-p "$THREAD" \
		--paired-end \
		--bam \
    --estimate-rspd \
		$EXP_DIR/4.BAM/"$SAMPLE_ID"."$GENOME_VERSION"_Aligned.toTranscriptome.out.bam \
		$REF_DIR/"$GENOME_VERSION" \
		$EXP_DIR/8.RSEM_geneCount/"$SAMPLE_ID"."$GENOME_VERSION" \
		2> $EXP_DIR/8.RSEM_geneCount/"$SAMPLE_ID"."$GENOME_VERSION".rsem.log

fi

rsem-plot-model \
	$EXP_DIR/8.RSEM_geneCount/"$SAMPLE_ID"."$GENOME_VERSION" \
	$EXP_DIR/8.RSEM_geneCount/"$SAMPLE_ID"."$GENOME_VERSION".FragLengthDist.pdf

else
	echo 'Step6: RSEM gene count (skip)'
fi

date

#############################################################################
################################### end #####################################
#############################################################################
