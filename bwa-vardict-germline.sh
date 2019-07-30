#!/usr/bin/env bash

echo '# Starting PYNEH NGS analysis pipeline version 20190228...'

# Check the required programs
echo '====='
echo "[Timestamp: `date`]"
echo '# Initializing...'
SCRIPT_DIR=`dirname "$(readlink -f "$0")"`
echo "    Script directory: $SCRIPT_DIR"
echo '== Programs =='
echo "    FastQC path: `which fastqc`"
echo "    BWA path: `which bwa`"
echo "    SPEEDSEQ path: `which speedseq`"
echo "    SAMTOOLS path: `which samtools`"
echo "    BCFTOOLS path: `which bcftools`"
echo "    CNVKIT path: `which cnvkit.py`"
echo "    VARDICT (core) path: `which VarDict`"
echo "    VARDICT (strand bias test) path: `which teststrandbias.R`"
echo "    R path: `which R`"
echo "    VARDICT (var2vcf) path: `which var2vcf_valid.pl`"
echo "    ANNOVAR path: `which table_annovar.pl`"
echo "    PYTHON3 path: `which python3`"
BUILD_BED_FILE=`readlink -e $SCRIPT_DIR/build_bed_file.py`
echo "    BUILD_BED_FILE path: $BUILD_BED_FILE"
IGV_PATH=`readlink -e ~/Programs/igv_utils/igv_plotter/lib/igv.jar`
echo "    IGV path: $IGV_PATH"
IGV_PLOTTER_PATH=`which igv_plotter`
echo "    IGV_PLOTTER path: $IGV_PLOTTER_PATH"
RUN_IGV_PLOTTER=`readlink -e $SCRIPT_DIR/RunIGVPlotter.sh`
echo "    RUN_IGV_PLOTTER path: $RUN_IGV_PLOTTER"
PDF_REPORT_TOOL=`readlink -e $SCRIPT_DIR/ngs_variant_report.py`
echo "    NGS_VARIANT_REPORT path: $PDF_REPORT_TOOL"
PDF_COVREPORT_TOOL=`readlink -e $SCRIPT_DIR/gene_panel_coverage_report.py`
echo "    NGS_COVERAGE_REPORT path: $PDF_COVREPORT_TOOL"
echo "    BEDTOOLS path: `which bedtools`"
echo '== Databases =='
HG19_PATH=`readlink -e ~/bundle/hg19/ucsc.hg19.fasta.gz`
if [[ -r $HG19_PATH ]]; then
	echo '    hg19... OK'
else
	echo '    hg19... Failure!'
	exit
fi

UNZIP_HG19_PATH=`readlink -e ~/bundle/hg19/ucsc.hg19.fasta`
if [[ -r $UNZIP_HG19_PATH ]]; then
	echo '    hg19 (unzip)... OK'
else
	echo '    hg19 (unzip)... Failure!'
	exit
fi

HUMANDB_PATH=`readlink -e ~/humandb`
if [[ -r $HUMANDB_PATH ]]; then
	echo '    humandb... OK'
else
	echo '    humandb... Failure!'
	exit
fi

CONTROL_BAM_PATH=`readlink -e ~/localdb/control.bam`
if [[ -r $CONTROL_BAM_PATH ]]; then
	echo '    control BAM file... OK'
else
	echo '    control BAM file... Failure!'
	exit
fi

SV_EXCLUDE_BED_PATH=`readlink -e ~/Programs/speedseq/annotations/ceph18.b37.exclude.2014-01-15.bed`
if [[ -r $SV_EXCLUDE_BED_PATH ]]; then
	echo '    SV exclusion BED file... OK'
else
	echo '    SV exclusion BED file... Failure!'
	exit
fi

# Show the help message if the required number of arguments is not found
if [[ $# -ne 8 ]]; then
	echo "    Usage: $0 [SAMPLE_NAME] [FASTQ1] [FASTQ2] [COVERAGE] [FLANK_BP] [GENE_LIST_TXT] [CNVKIT_REF] [NUM_THREADS]"
	echo $#
	exit
fi

# Check the input files
echo '# Validating inputs...'
echo "    Sample name: $1"

# Determine the sample type
SAMPLE_TYPE='undefined'
if [[ $1 == *"WES"* ]]; then
	SAMPLE_TYPE='WES'
	echo "    ..Sample type => WES"
fi
if [[ $1 == *"WGS"* ]]; then
	SAMPLE_TYPE='WGS'
	echo "    ..Sample type => WGS"
fi
if [[ $SAMPLE_TYPE == 'undefined' ]]; then
	echo "    ..Sample type could not be inferred from $1!"
        exit
fi

# Determine if CNVkit analysis should be performed
DO_CNVKIT='undefined'
if [[ $7 == 'NONE' ]]; then
        DO_CNVKIT=false
        echo "    ..CNVkit analysis will NOT be performed."
fi
if [[ -r $7 ]]; then
        DO_CNVKIT=true
	CNVREF_PATH=$7
        echo "    ..CNVkit analysis will be performed using reference $CNVREF_PATH."
fi
if [[ $SAMPLE_TYPE == 'undefined' ]]; then
	echo "    ..CNVkit analysis status could not be inferred from $7!"
        exit
fi

echo "    FASTQ 1: $2"
if [[ -r $2 ]]; then
	echo "    ..File check... OK"
	FASTQ1=$2
else
	echo "    ..File check... Failure!"
	exit
fi

echo "    FASTQ 2: $3"
if [[ -r $3 ]]; then
	echo "    ..File check... OK"
	FASTQ2=$3
else
	echo "    ..File check... Failure!"
	exit
fi

DIR_A=`dirname $2`
DIR_B=`dirname $3`
if [[ "$DIR_A" != "$DIR_B" ]]; then
	echo "    FASTQ files in different directories!"
	exit
else
	echo "# Setting output directory to: $DIR_A"
	OUTPUT_DIR=$DIR_A
fi

COVERAGE=$4
echo "    Minimum coverage set to: $COVERAGE x"

FLANK_BP=$5
echo "    CDS flanking set to: $FLANK_BP bp"

echo "    Gene list text file: $6"
if [[ -r $6 ]]; then
	echo "    ..File check... OK"
	GENE_LIST_TXT=$6
else
	echo "    ..File check... Failure!"
	exit
fi

# Set the number of threads
NUM_THREADS=$8
echo "Analysis will start using $NUM_THREADS threads"

# Perform mapping
BAM_PATH="$OUTPUT_DIR/$1.bam"
BAI_PATH="$OUTPUT_DIR/$1.bam.bai"
SPLITTERS_PATH="$OUTPUT_DIR/$1.splitters.bam"
DISCORDANTS_PATH="$OUTPUT_DIR/$1.discordants.bam"
if [ -r $BAM_PATH ] && [ -r $BAI_PATH ] && [ -r $SPLITTERS_PATH ] && [ -r $DISCORDANTS_PATH ]; then
	echo "[Timestamp: `date`]"
	echo '!! Files found. Skipping Step 1. !!'
else
	echo '====='
	echo "[Timestamp: `date`]"
	echo '# Step 1: Mapping, deduplication and sorting...'
	echo "# Writing output to $BAM_PATH"
	speedseq align -R "@RG\tID:$1\tSM:$1\tLB:$1\tPL:ILMN_LIKE" -t $NUM_THREADS -o $OUTPUT_DIR/$1 $HG19_PATH $FASTQ1 $FASTQ2
fi

# Generate BED file
BED_100_PATH="$OUTPUT_DIR/$1.flankplus100.bed"
PLUS_100_FLANK=$((FLANK_BP+100))
BED_PATH="$OUTPUT_DIR/$1.bed"
if [ -r $BED_PATH ]; then
	echo "[Timestamp: `date`]"
	echo '!! File found. Skipping Step 2. !!'
else
	echo '====='
	echo "[Timestamp: `date`]"
	echo '# Step 2: Generate BED files for variant calling...'
	echo "# Writing output to $BED_100_PATH and $BED_PATH"
	GENE_STRING=`grep -v '^#'  $GENE_LIST_TXT | grep -v '^$' | tr '\r\n' '\n' | tr '\n' ' '`
	python3 $BUILD_BED_FILE -f $PLUS_100_FLANK $GENE_STRING > $BED_100_PATH
	python3 $BUILD_BED_FILE -f $FLANK_BP $GENE_STRING > $BED_PATH
fi

# Perform small indel and SNV calling
PLUS_100_VCF_PATH="$OUTPUT_DIR/$1.rawplus100.vcf"
RAW_VCF_PATH="$OUTPUT_DIR/$1.raw.vcf"
if [ -r $RAW_VCF_PATH ]; then
	echo "[Timestamp: `date`]"
	echo '!! File found. Skipping Step 3. !!'
else
	echo '====='
	echo "[Timestamp: `date`]"
	echo '# Step 3: Small variant calling...'
	echo "# Writing output to $RAW_VCF_PATH"
	VarDict -G $UNZIP_HG19_PATH -f 0.05 -I 1000 -L 1001 -k 1 -N $1 -th $NUM_THREADS --dedup -b $BAM_PATH -z -c 1 -S 2 -E 3 -g 4 $BED_100_PATH | teststrandbias.R | var2vcf_valid.pl -N $1 -E -f 0.05 > $PLUS_100_VCF_PATH
	echo '# Intersecting additionally flanked output with BED file...'
	bedtools intersect -header -a $PLUS_100_VCF_PATH -b $BED_PATH > $RAW_VCF_PATH
fi

# Perform small indel and SNV annotation
ANNO_VCF_PREFIX="$OUTPUT_DIR/$1.annovar"
ANNO_VCF_PATH="$OUTPUT_DIR/$1.annovar.hg19_multianno.vcf"
if [ -r $ANNO_VCF_PATH ]; then
	echo "[Timestamp: `date`]"
	echo '!! File found. Skipping Step 4. !!'
else
	echo '====='
	echo "[Timestamp: `date`]"
	echo '# Step 4: Variant annotation...'
	echo "# Writing output to $ANNO_VCF_PATH"
	table_annovar.pl $RAW_VCF_PATH $HUMANDB_PATH -buildver hg19 -out $ANNO_VCF_PREFIX -remove -protocol refGene,ensGene,1000g2015aug_all,1000g2015aug_eas,gnomad_genome,exac03,gnomad_exome,dbnsfp35c,dbscsnv11,intervar_20180118,clinvar_20180603,avsnp150 -operation g,g,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput
fi

# Perform SV calling
SV_VCF_PATH="$OUTPUT_DIR/$1.sv.vcf.gz"

if [ -r $SV_VCF_PATH ]; then
	echo "[Timestamp: `date`]"
	echo '!! File found. Skipping Step 5. !!'
else
	echo '====='
	echo "[Timestamp: `date`]"
	echo '# Step 5: Structural variant calling...'
	echo "# Writing output to $SV_VCF_PATH"
        if [[ $SAMPLE_TYPE == 'WES' ]]; then
		# LUMPY analysis only
		speedseq sv -B $BAM_PATH -S $SPLITTERS_PATH -D $DISCORDANTS_PATH -R $HG19_PATH -o $OUTPUT_DIR/$1 -x $SV_EXCLUDE_BED_PATH -t $NUM_THREADS -P
	fi
	if [[ $SAMPLE_TYPE == 'WGS' ]]; then
		# LUMPY + CNVnator for WGS
		speedseq sv -B $BAM_PATH -S $SPLITTERS_PATH -D $DISCORDANTS_PATH -R $HG19_PATH -o $OUTPUT_DIR/$1 -x $SV_EXCLUDE_BED_PATH -t 8 -d -P
		echo "# Intersecting region of ROI bed file with read depth SV calls..."
		RD_BED_PATH="$OUTPUT_DIR/$1.sv.$1.bam.readdepth.bed"
		bedtools intersect -a $BED_PATH -b $RD_BED_PATH -wb > $OUTPUT_DIR/$1.ROI.SV_report.txt
	fi
fi

# Perform IGV plotting

echo '====='
echo "[Timestamp: `date`]"
echo '# Step 6: Automatic IGV plotting...'
echo "# Writing output to $OUTPUT_DIR/snapshots"

if [ -r "$OUTPUT_DIR/snapshots" ]; then
	echo "IGV snapshot directory found."
	echo '!! File found. Skipping Step 6. !!'
else
	echo "Creating IGV snapshot directory..."
	mkdir $OUTPUT_DIR/snapshots
	cd "$OUTPUT_DIR/snapshots"
	echo "Current working directory: `pwd`"
	$RUN_IGV_PLOTTER $1 $BAM_PATH $BED_PATH $ANNO_VCF_PATH $CONTROL_BAM_PATH $SV_VCF_PATH $IGV_PATH $IGV_PLOTTER_PATH
	cd $SCRIPT_DIR
	echo "Current working directory: `pwd`"
fi

# Compile PDF variant report

echo '====='
echo "[Timestamp: `date`]"
echo '# Step 7: Compile PDF variant report'
if [ -r "$OUTPUT_DIR/$1.pdf" ]; then
	echo "Main PDF report found."
	echo '!! File found. Skipping Step 7. !!'
else
	echo "# Output will be written to $OUTPUT_DIR/$1.pdf"
	python3 $PDF_REPORT_TOOL $ANNO_VCF_PATH
	if [ -r "$ANNO_VCF_PATH.pdf" ]; then
		echo 'Renaming output...'
		mv "$ANNO_VCF_PATH.pdf" "$OUTPUT_DIR/$1.pdf"
	else
		echo 'Output file not found!'
		exit
	fi
fi


# Compile PDF coverage report

echo '====='
echo "[Timestamp: `date`]"
echo '# Step 8: Compile PDF coverage report'
if [ -r "$OUTPUT_DIR/$1.coverage.pdf" ]; then
	echo "PDF coverage report found."
	echo '!! File found. Skipping Step 8. !!'
else
	echo "# Output will be written to $OUTPUT_DIR/$1.coverage.pdf"
	python3 $PDF_COVREPORT_TOOL $BAM_PATH $COVERAGE $FLANK_BP $GENE_LIST_TXT
	if [ -r "$BAM_PATH.coverage.pdf" ]; then
		echo 'Renaming output...'
		mv "$BAM_PATH.coverage.pdf" "$OUTPUT_DIR/$1.coverage.pdf"
	else
		echo 'Output file not found!'
		exit
	fi
fi


# Optionally, perform the force-call of SNPs

echo '====='
echo "[Timestamp: `date`]"
echo '# Step 9: Perform force-call of variants'
if [ -r "$OUTPUT_DIR/FORCECALL.snp" ]; then
	echo 'Force-call target file found.'
	echo 'Now building force-call region list...'
	TARGET_STRING=`grep -v '^#' "$OUTPUT_DIR/FORCECALL.snp" | grep -v '^$' | tr '\r\n' '\n' | tr '\n' ' '`
	python3 $BUILD_BED_FILE -f $FLANK_BP $TARGET_STRING > "$OUTPUT_DIR/FORCECALL.regions"
	echo 'Now performing force-calling...'
	bcftools mpileup -f $UNZIP_HG19_PATH -R $OUTPUT_DIR/FORCECALL.regions $BAM_PATH | bcftools call -m > $OUTPUT_DIR/$1.FORCECALL.vcf
	echo 'Now performing additional IGV plotting...'
	cd "$OUTPUT_DIR/snapshots"
	echo "Current working directory: `pwd`"
	if [[ $SAMPLE_TYPE == 'WGS' ]]; then
		$RUN_IGV_PLOTTER $1 $BAM_PATH $BED_PATH $OUTPUT_DIR/$1.FORCECALL.vcf $CONTROL_BAM_PATH $RD_BED_PATH $IGV_PATH $IGV_PLOTTER_PATH
	fi
	if [[ $SAMPLE_TYPE == 'WES' ]]; then
		$RUN_IGV_PLOTTER $1 $BAM_PATH $BED_PATH $OUTPUT_DIR/$1.FORCECALL.vcf $CONTROL_BAM_PATH $SV_VCF_PATH $IGV_PATH $IGV_PLOTTER_PATH
	fi
	cd $SCRIPT_DIR
	echo "Current working directory: `pwd`"
	echo 'Now performing additional reporting...'
	echo "# Output will be written to $OUTPUT_DIR/$1.FORCECALL.pdf"
	python3 $PDF_REPORT_TOOL $OUTPUT_DIR/$1.FORCECALL.vcf $OUTPUT_DIR/FORCECALL.regions
	if [ -r "$OUTPUT_DIR/$1.FORCECALL.vcf.pdf" ]; then
		echo 'Renaming output...'
		mv "$OUTPUT_DIR/$1.FORCECALL.vcf.pdf" "$OUTPUT_DIR/$1.FORCECALL.pdf"
	else
		echo 'Output file not found!'
		exit
	fi
else
	echo "Force-call target file ($OUTPUT_DIR/FORCECALL.snp) not found. Skipping force-calling.."
fi

# Optionally, perform CNVkit analysis
echo '====='
echo "[Timestamp: `date`]"
echo '# Step 10: Perform CNVkit analysis'
if [ $DO_CNVKIT == true ] && [ -r "$CNVREF_PATH" ]; then
	echo "CNVkit reference file found."
	cnvkit.py batch $BAM_PATH -r $CNVREF_PATH --output-dir $OUTPUT_DIR/cnvkit -p $NUM_THREADS
	echo "# Intersecting CNVkit output with custom BED file..."
	grep -v 'chromosome' $OUTPUT_DIR/cnvkit/$1.cnr > $OUTPUT_DIR/cnvkit/$1.noheader.cnr
	bedtools intersect -wb -a $BED_PATH -b $OUTPUT_DIR/cnvkit/$1.noheader.cnr > $OUTPUT_DIR/cnvkit.ROI.txt
else
	echo "CNVkit reference file not found. Skipping CNVkit analysis.."
fi

# Generate quality statistics

echo '====='
echo "[Timestamp: `date`]"
echo '# Step 11: Generate FastQC reports'
echo "# Output will be written to $OUTPUT_DIR/fastqc"
if [ -d "$OUTPUT_DIR/fastqc" ]; then
	echo '!! FastQC output directory found. Skipping Step 11. !!'
	echo "Note: Remove $OUTPUT_DIR/fastqc in order to trigger the FastQC step."
else
	cd $OUTPUT_DIR
	echo "Current working directory: `pwd`"
	fastqc -t $NUM_THREADS $FASTQ1 $FASTQ2
	mkdir fastqc
	mv *.html *.zip ./fastqc
	cd $SCRIPT_DIR
	echo "Current working directory: `pwd`"
fi

# Generate ROI BAM file for storage

echo '====='
echo "[Timestamp: `date`]"
echo '# Step 12: Generate ROI BAM files for storage'
echo "# Output will be written to $OUTPUT_DIR/$1.ROI.bam"
if [ -r "$OUTPUT_DIR/$1.ROI.bam" ]; then
	echo 'ROI BAM file found. Step 10 will be skipped.'
	echo "Note: Remove $OUTPUT_DIR/$1.ROI.bam in order to trigger the ROI BAM generation step."
else
	bedtools intersect -a $BAM_PATH -b $BED_PATH > $OUTPUT_DIR/$1.ROI.bam
	samtools index $OUTPUT_DIR/$1.ROI.bam
fi

echo '====='
echo 'DONE!'
echo "[Timestamp: `date`]"

