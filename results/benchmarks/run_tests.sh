#!/bin/bash
# IMPORTANT: must be run from output/mm9 or similar

DATA_DIR="${PWD}/../../data"
ROOT_DIR=$PWD
EXEC_CMD="${PWD}/../../triplexator/bin/triplexator"
PALINDROM=false
MYERS=false
BRUTE=false
QGRAM=false

l=25
eLow=5
eHigh=20
cLow=1
cHigh=3

while getopts ":mbqpl:c:t:" opt; do
    case ${opt} in
        b ) BRUTE=true
        ;;
        m ) MYERS=true
        ;;
        q ) QGRAM=true
        ;;
        l ) l=$OPTARG
        ;;
        c ) chr=$OPTARG
        ;;
        \? ) echo "Usage: run_tests.sh [-m] [-i] [-l min-length]"
        ;;
    esac
done

DNA_FILE="${DATA_DIR}/dna/hg38/genome_hg38_chr1.fa"
#RNA_FILE="${DATA_DIR}/rna/D0_D4_up/ENSG00000_1_merged.fa"
RNA_FILE="${DATA_DIR}/rna/D0_D4_up/ENSG_merged_complete.fa"
#OUTPUT_DIR="D0_D4_up_ENSG00000_1__hg38_chr1"
OUTPUT_DIR="D0_D4_up_ENSG00000_complete__hg38_chr1"

MYERS_TMP="${ROOT_DIR}/myers.tmp"
BRUTE_TMP="${ROOT_DIR}/brute.tmp"
OUTPUT_PREFIX=""

for ((e=$eLow;e<=$eHigh;e+=5)); do
    for ((c=$cLow;c<=$cHigh;c+=1)); do
        if [ "$MYERS" = true ]; then
            cd ${ROOT_DIR}/${OUTPUT_DIR}/bit_parallel/
            OUTPUT="${OUTPUT_PREFIX}l${l}_c${c}_e${e}"

            echo "Running l: ${l}, c: ${c}, e: ${e}"
            ulimit -Sv 25000000 
            $EXEC_CMD -ss $RNA_FILE -ds $DNA_FILE -o ${OUTPUT}.tpx --bit-parallel -fm 1 -v -l $l -c $c -e $e &> ${OUTPUT}.dbg
            #sort ${OUTPUT}.tpx > ${MYERS_TMP}
        fi

        if [ "$BRUTE" = true ]; then
#TODO revert to normal / not palindrom
            cd ${ROOT_DIR}/palindrom
            if [[ "$TRUNCATED" = 0 ]]; then
                DNA=$DNA_FILE
                RNA=$RNA_FILE
                OUTPUT="chr${chr}_brute_l${l}_c${c}_e${e}"
            elif [[ "$TRUNCATED" = 1 ]]; then
                DNA=$DNA_FILE
                RNA=$RNA_TRUNCATED_FILE
                OUTPUT="chr${chr}_brute_l${l}_c${c}_e${e}_truncated_fulldna"
            elif [[ "$TRUNCATED" = 2 ]]; then
                DNA=$DNA_TRUNCATED_FILE
                RNA=$RNA_TRUNCATED_FILE
                OUTPUT="chr${chr}_brute_l${l}_c${c}_e${e}_truncated"
            fi

            echo "Running BRUTE l: ${l}, c: ${c}, e: ${e}"
#    echo "$EXEC_CMD -ss $RNA -ds $DNA -o ${OUTPUT}.tpx -v -l $l -c $c -e $e &> ${OUTPUT}.dbg"
            $EXEC_CMD -ss $DNA -ds $DNA -o ${OUTPUT}.tpx -v -l $l -c $c -e $e &> ${OUTPUT}.dbg
            sort ${OUTPUT}.tpx > ${BRUTE_TMP}
        fi

        if [ "$QGRAM" = true ]; then
            cd ${ROOT_DIR}/qgram
            if [[ "$TRUNCATED" = 0 ]]; then
                DNA=$DNA_FILE
                RNA=$RNA_FILE
                OUTPUT="chr${chr}_qgram_l${l}_c${c}_e${e}_q2"
            fi

            echo "Running l: ${l}, c: ${c}, e: ${e} q?"
            /home/zbarni/code/triplexator-project/triplexator_original/bin/triplexator -ss $RNA -ds $DNA -o ${OUTPUT}.tpx -v -t 1 -fm 1 -l $l -c $c -e $e &> ${OUTPUT}.dbg
        fi
    done
done

#if [ ! -f ${MYERS_TMP} ]; then
#    echo "Myers sorted output not found. Can't compare."
#elif [ ! -f ${BRUTE_TMP} ]; then
#    echo "Brute sorted output not found. Can't compare."
#else
#    echo "Comparing myers and brute-force... wait for failed?"
#    cmp --silent ${MYERS_TMP} ${BRUTE_TMP} || echo "Test #${TEST} failed."
#fi
#rm -f ${MYERS_TMP} ${BRUTE_TMP}
