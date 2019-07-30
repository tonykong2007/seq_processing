#!/bin/bash

if [[ $# -ne 8 ]] ; then
    echo 'igv_plotter wrapper script by Tom C.C. Ho (c) 2017-2019'
    echo '-----'
    echo "Usage: RunIGVPlotter.sh [SAMPLE_NAME] [SAMPLE_BAM_PATH] [ROI_BED_PATH] [SAMPLE_VCF_PATH] [CONTROL_BAM_PATH] [SAMPLE_SV_VCF_PATH] [IGV_PATH] [IGV_PLOTTER_PATH]"
    echo ''
    exit 0
fi

# essential file paths
PATH_TO_TEMP_FILE="./temp.loci"
PATH_TO_SAMPLE=$2
PATH_TO_ROI=$3
PATH_TO_VCF=$4
PATH_TO_CONTROL=$5
PATH_TO_SV_VCF=$6
PATH_TO_IGV=$7
PATH_TO_IGVPLOTTER=$8

# essential parameters

# generate target regions for plotting
grep "^[^#]" $PATH_TO_VCF | cut -f1,2 --output-delimiter ':' > $PATH_TO_TEMP_FILE

# perform automated IGV plotting
echo "Now plotting sample ${1}..."
$PATH_TO_IGVPLOTTER --igv-jar-path $PATH_TO_IGV -p SAM.SHOW_CENTER_LINE=TRUE -p SAM.ALLELE_THRESHOLD=0.005 -p SAM.ALLELE_USE_QUALITY=FALSE -p SAM.SHOW_SOFT_CLIPPED=TRUE -p SAM.COLOR_BY=READ_STRAND -m 16G --squish -v -L $PATH_TO_TEMP_FILE -o $1 $PATH_TO_SAMPLE $PATH_TO_CONTROL $PATH_TO_VCF $PATH_TO_SV_VCF $PATH_TO_ROI
echo "Now removing temp file..."
rm $PATH_TO_TEMP_FILE
