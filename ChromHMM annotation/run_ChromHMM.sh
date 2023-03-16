#!/bin/bash

set -e

COMMAND=$1

##### Usage
usage(){
  echo -e "Usage: bash ChromHMM.sh [BinarizeBam|LearnModel]
                  Example: bash ChromHMM.sh BinarizeBam INPUT_DIR SAMPLE_TAB OUTPUT_DIR
		           bash ChromHMM.sh BinarizeBed INPUT_DIR SAMPLE_TAB OUTPUT_DIR
                           bash ChromHMM.sh LearnModel BinarizeDIR OUTPUT_DIR NUM_MODEL
               " 1>&2
}
exit_abnormal(){
  usage
  exit 1
}

if [[ $# != 4 ]]
then
    echo -e "Error: One argument (either \"BinarizeBam\" or \"LearnModel\") is requred!!!"
    exit_abnormal
    exit 1
fi


if [[ ${COMMAND} == "BinarizeBam" ]]
then
    INPUT_DIR=$2
    SAMPLE_TAB=$3
    OUTPUT_DIR=$4
    java -Xmx12G -XX:ParallelGCThreads=24 -jar ChromHMM/ChromHMM.jar ${COMMAND} -p 0.000001 -f 4 -gzip -p 24 refGenome/Chromosome_Lengths.txt ${INPUT_DIR} ${SAMPLE_TAB} ${OUTPUT_DIR}/
elif [[ ${COMMAND} == "BinarizeBed" ]]
then
    INPUT_DIR=$2
    SAMPLE_TAB=$3
    OUTPUT_DIR=$4
    java -Xmx12G -XX:ParallelGCThreads=24 -jar ChromHMM/ChromHMM.jar ${COMMAND} -p 24 refGenome/Chromosome_Lengths.txt ${INPUT_DIR} ${SAMPLE_TAB} ${OUTPUT_DIR}/
elif [[ ${COMMAND} == "LearnModel" ]]
then
    BinarizeDIR=$2
    OUTPUT_DIR=$3
    NUM_MODEL=$4
    java -Xmx24G -XX:ParallelGCThreads=24 -jar ChromHMM/ChromHMM.jar ${COMMAND} -p 24 -l refGenome/Chromosome_Lengths.txt ${BinarizeDIR} ${OUTPUT_DIR} ${NUM_MODEL} bosTau9
fi
