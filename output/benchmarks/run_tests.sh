#!/bin/bash
# IMPORTANT: must be run from output/mm9 or similar

l=15
chr=1
DATA_DIR="${PWD}/../../data"
ROOT_DIR=$PWD/..
EXEC_CMD="${PWD}/../../triplexator/bin/triplexator"
MYERS=false
BRUTE=false
TRUNCATED=0
eLow=5
eHigh=5
cLow=2
cHigh=2

while getopts ":mbl:c:t:" opt; do
    case ${opt} in
        t ) TRUNCATED=$OPTARG
        ;;
        b ) BRUTE=true
        ;;
        m ) MYERS=true
        ;;
        l ) l=$OPTARG
        ;;
        c ) chr=$OPTARG
        ;;
        \? ) echo "Usage: run_tests.sh [-m] [-i] [-l min-length]"
        ;;
    esac
done

DNA_FILE="${DATA_DIR}/dna/mm9/mm9.chr${chr}.fa"
DNA_TRUNCATED_FILE="${DATA_DIR}/dna/mm9/mm9.chr${chr}.oneline.truncated.fa"
RNA_FILE="${DATA_DIR}/rna/CDP_merged.fa"
RNA_TRUNCATED_FILE="${DATA_DIR}/rna/CDP_truncated.fa"

MYERS_TMP="${ROOT_DIR}/myers.tmp"
BRUTE_TMP="${ROOT_DIR}/brute.tmp"

for ((e=$eLow;e<=$eHigh;e+=5)); do
    for ((c=$cLow;c<=$cHigh;c+=1)); do

# switch directory according to options
        if [ "$MYERS" = true ]; then
            cd ${ROOT_DIR}/myers
            if [[ "$TRUNCATED" = 0 ]]; then
                DNA=$DNA_FILE
                RNA=$RNA_FILE
                OUTPUT="chr${chr}_myers_l${l}_c${c}_e${e}"
            elif [[ "$TRUNCATED" = 1 ]]; then
                DNA=$DNA_FILE
                RNA=$RNA_TRUNCATED_FILE
                OUTPUT="chr${chr}_myers_l${l}_c${c}_e${e}_truncated_fulldna"
            elif [[ "$TRUNCATED" = 2 ]]; then
                DNA=$DNA_TRUNCATED_FILE
                RNA=$RNA_TRUNCATED_FILE
                OUTPUT="chr${chr}_myers_l${l}_c${c}_e${e}_truncated"
            fi

            echo "Running l: ${l}, c: ${c}, e: ${e}"
#    echo "$EXEC_CMD -ss $RNA -ds $DNA -o ${OUTPUT}.tpx -i -v -l $l -c $c -e $e &> ${OUTPUT}.dbg"
            $EXEC_CMD -ss $RNA -ds $DNA -o ${OUTPUT}.tpx -i -fm 1 -v -l $l -c $c -e $e &> ${OUTPUT}.dbg

            sort ${OUTPUT}.tpx > ${MYERS_TMP} 
            #rm ${OUTPUT}.debug
        fi

        if [ "$BRUTE" = true ]; then
            cd "${ROOT_DIR}/mm9/original/"
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

            echo "Running l: ${l}, c: ${c}, e: ${e}"
#    echo "$EXEC_CMD -ss $RNA -ds $DNA -o ${OUTPUT}.tpx -v -l $l -c $c -e $e &> ${OUTPUT}.dbg"
            $EXEC_CMD -ss $RNA -ds $DNA -o ${OUTPUT}.tpx -v -l $l -c $c -e $e &> ${OUTPUT}.dbg
            sort ${OUTPUT}.tpx > ${BRUTE_TMP}
        fi
    done
done

if [ ! -f ${MYERS_TMP} ]; then
    echo "Myers sorted output not found. Can't compare."
elif [ ! -f ${BRUTE_TMP} ]; then
    echo "Brute sorted output not found. Can't compare."
else
    echo "Comparing myers and brute-force... wait for failed?"
    cmp --silent ${MYERS_TMP} ${BRUTE_TMP} || echo "Test #${TEST} failed."
fi
#rm -f ${MYERS_TMP} ${BRUTE_TMP}
