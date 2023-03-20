#!/bin/bash -l
set -e
# =====================================================================
# Generate the signal track files in the ROADMAP/ENCODE processing pipeline (for imputation)
#   1. Fold change over control at base resolution
#   2. -log10p.val by Poisson model at base resolution
# 
#   Adapted from Anshul Kundaje. See "https://sites.google.com/site/anshulkundaje/projects/epigenomeroadmap#TOC-Uniformly-processed-and-normalized-genome-wide-signal-coverage-tracks" for more in-depth explanation of this pipeline.
# NOTE: Compares a single chip-seq dataset against a single control.
# =====================================================================
library=$1
echo -e "Processing ${library}...\n\n"

# Arguments / Software:
GENOMESIZE="2715853792"
GENOME='/group/zhougrp3/dguan/BovineFAANG/refGenome/Chromosome_Lengths.txt' 
# Check if MACS2.1 is in path: 
# (load with reuse .macs2-2.1.1.20160309-python-2.7.1-sqlite3-rtrees)
source ~/.bashrc
conda activate macs2
if [[ -z $(which macs2) ]]; then echo 'ERROR: MACS executable not in $PATH'; exit 1; fi
# Alternative MACS2.0 path: MACSPATH='/broad/compbio/anshul/projects/encode/preprocessing/peakcalling/macs/scripts/macs2'

# ========================
# Data and Pre-processing:
# ========================
# Check file existence:
if [ ! -f Aligned_Reads/${library}.bam ]; then echo 'ERROR: chip bam file not exists!'; exit 1; fi
control="Input_${library#*_}"
if [[ ! ${library} =~ ATAC* ]]
then
    if [ ! -f Aligned_Reads/${control}.bam ]; then echo 'ERROR: control bam file not exists!'; exit 1; fi
fi
if [ ! -f Metrics/${library}.spp_stats.txt ]; then echo 'ERROR: spp file not exists!'; exit 1; fi
if [ ! -f Aligned_Reads/${library}.tagAlign.gz ]; then echo 'ERROR: chip tag file not exists!'; exit 1; fi
if [[ ! ${library} =~ ATAC* ]]
then
    if [ ! -f Aligned_Reads/${control}.tagAlign.gz ]; then echo 'ERROR: control tag not exists!'; exit 1; fi
fi

# Output file name
PKPREF="Bedgraphs/${library}"
if [[ -s "${PKPREF}.pval.signal.bedgraph" ]]
then
    echo "DONE: Skipping ${PKPREF} as it is already done."
    exit 1
fi

# ==================================
# Script for signal track generation
# (Previously submitted separately).
# ==================================
tmpdir="Temp"
mkdir -p ${tmpdir}

# --------------------------
# Determine if reads are paired
# or single
# --------------------------
#nreads=`ls -1 Raw_Reads/${library}*.fq.gz | wc -l`
#if [[ ${nreads} -eq 1 ]]
#then
#    file_format="-f BAM"
#elif [[ ${nreads} -eq 2 ]]
#then
#    file_format="-f BAMPE"
#else
#    echo "${library} is either single or peaired reads, please manually check..."
#fi

nreads=$(cat Seq_file_type.txt | awk -v i=${library} '{if ($1 == i) print $2}')
if [[ ${nreads} == "Paired" ]]
then
    file_format="-f BAMPE"
elif [[ ${nreads} == "Single" ]]
then
    file_format="-f BAM"
else
    echo "${library} is either single or peaired reads, please manually check..."
fi
echo -e "ChIP library is: ${nreads}...\n\n"

# --------------------------
# Get read depth ratio and
# create tmp copy of control
# --------------------------
chiptag="Aligned_Reads/${library}.tagAlign.gz"
inputtag="Aligned_Reads/${control}.tagAlign.gz"
chip_depth=$(zcat ${chiptag} | wc -l)
echo -e "Number of ChIP reads: ${chip_depth}\n"
if [[ ! ${library} =~ ATAC* ]]
then
    input_depth=$(zcat ${inputtag} | wc -l)
else
   input_depth=0
   echo -e "Number of Input reads: ${chip_depth}\n"
fi

if [[ ${chip_depth} -gt ${input_depth} || ! ${library} =~ ATAC* ]]
then
    scaling=" --to-large "
else
    scaling=""
fi
chipReads=$(zcat ${chiptag} | wc -l | awk '{printf "%f", $1/1000000}') 
if [[ ! ${library} =~ ATAC* ]]
then
    controlReads=$(zcat ${inputtag} | wc -l | awk '{printf "%f", $1/1000000}' )
    #echo -e "Number of ChIP reads: ${chipReads} (x 10^6)\nNumber of Input reads: ${controlReads} (x 10^6)\n\n\n"
    sval=$(echo "${chipReads} ${controlReads}" | awk '$1>$2{printf "%f",$2} $1<=$2{printf "%f",$1}' )
else
	sval=1.0
fi
echo -e "Scaling factor for treatment and control track: $sval ... \n\n\n"

# ===================================================================
# Create the pileup and control lambda bedgraph tracks using MACS2.1:
# ===================================================================
if [ ! -f "${PKPREF}_control_lambda.bdg" ] || [ ! -f "${PKPREF}_treat_pileup.bdg" ]
then
    if [[ ${library} == H3K27me3* || ${library} == H3K36me3* || ${library} == H3K9me3* ]]
    then
        libtype="--broad"
        peak_file="${PKPREF}_peaks.broadPeak"
        combchip="-t Aligned_Reads/${library}.bam"
        combcontrol="-c Aligned_Reads/${control}.bam --fix-bimodal"
        pvalue="-q 0.05"
    elif [[ ${library} == H3K4me3* || ${library} == H3K27ac* || ${library} == H3K4me1* || ${library} == CTCF*  ]]
    then
        libtype=""
        peak_file="${PKPREF}_peaks.narrowPeak"
        combchip="-t Aligned_Reads/${library}.bam"
        combcontrol="-c Aligned_Reads/${control}.bam --fix-bimodal "
        pvalue="-q 0.01"
    elif [[ ${library} == ATAC* ]]
    then
        libtype=""
        peak_file="${PKPREF}_peaks.narrowPeak"
        combchip="-t Aligned_Reads/${library}.bam --shift -100 --nomodel"
        combcontrol=""
        pvalue=""
    else
        echo -e "${library} is not chipseq, please manually check..."
    fi
    macs2 callpeak ${combchip} ${combcontrol} ${file_format} -n ${PKPREF} -g ${GENOMESIZE} ${pvalue} --extsize 200 -B --keep-dup all --SPMR ${libtype} ${scaling}
    gzip -f -c ${peak_file} > ${peak_file}.gz
    rm -f ${peak_file}
    rm -f ${PKPREF}_peaks.xls
    rm -f ${PKPREF}_summits.bed
else
    echo -e "${library} lambda bedgraph files are already there!!!\n"
fi

# ===================
# Generate bedgraphs:
# ===================
# FoldChange: 
FCPREF=${PKPREF}.fc.signal
nPeaks=`zcat ${peak_file}.gz | wc -l`
#if [[ ${nPeaks} -gt 8000 ]]
#then
  if [[ ! -e "${FCPREF}.bedgraph.gz" ]]
  then
    echo -e "Generating FoldChange bedgraph...\n"
    macs2 bdgcmp -t ${PKPREF}_treat_pileup.bdg -c ${PKPREF}_control_lambda.bdg -o ${PKPREF}_FE.bdg -m FE
    slopBed -i ${PKPREF}_FE.bdg -g ${GENOME} -b 0 | bedClip stdin ${GENOME} ${FCPREF}.bedgraph
    echo -e "Sort and create FC bigWig file"
    sort -k1,1 -k2,2n ${FCPREF}.bedgraph > ${FCPREF}.bedgraph.tmp
    # bedGraphToBigWig ${FCPREF}.bedgraph.tmp ${GENOME} ${LVPREF}.bigwig
    gzip -c ${FCPREF}.bedgraph.tmp > ${FCPREF}.bedgraph.gz
    rm -f ${FCPREF}.bedgraph.tmp
    rm -f ${FCPREF}.bedgraph
    rm -f ${PKPREF}_FE.bdg
  fi

  # -log10pval: 
  LVPREF=${PKPREF}.pval.signal
  if [[ ! -e "${LVPREF}.bedgraph.gz" ]]
  then
    echo -e "Generating log10pval bedgraph...\n"
    macs2 bdgcmp -t ${PKPREF}_treat_pileup.bdg -c ${PKPREF}_control_lambda.bdg -o ${PKPREF}_ppois.bdg -m ppois -S ${sval}
    slopBed -i ${PKPREF}_ppois.bdg -g ${GENOME} -b 0 | bedClip stdin ${GENOME} ${LVPREF}.bedgraph

    echo -e "Sort and create log10pval bigWig file...\n"
    sort -k1,1 -k2,2n ${LVPREF}.bedgraph > ${LVPREF}.bedgraph.tmp
    # bedGraphToBigWig ${LVPREF}.bedgraph.tmp ${GENOME} ${LVPREF}.bigwig
    gzip -c -f ${LVPREF}.bedgraph.tmp > ${LVPREF}.bedgraph.gz
    rm -f ${LVPREF}.bedgraph.tmp
    rm -f ${LVPREF}.bedgraph
    rm -f ${PKPREF}_ppois.bdg
  fi
#fi
if [[ -f ${LVPREF}.bedgraph.gz ]] && [[ -f ${FCPREF}.bedgraph.gz ]]
then
  for file in ${PKPREF}_control_lambda.bdg ${PKPREF}_treat_pileup.bdg;
  do
    gzip -c -f ${file} > ${file}.gz
    rm -f ${file}
  done
fi


if [[ -f ${PKPREF}_peaks.gappedPeak ]]
then
    rm -f ${PKPREF}_peaks.gappedPeak
fi

echo "Done!"
