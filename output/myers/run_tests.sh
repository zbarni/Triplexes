#!/bin/bash
# IMPORTANT: must be run from output/mm9 or similar

mem=false   #memory profiling
l=25
chr=1
DATA_DIR="${PWD}/../../data"
ROOT_DIR=$PWD
EXEC_CMD="${PWD}/../../triplexator/bin/triplexator"
LOCAL=true
DEBUG=false
MYERS=false
eLow=15
eHigh=15
cLow=1
cHigh=1

while getopts ":mdoil:c:" opt; do
    case ${opt} in
        m ) MYERS=true
        ;;
        d ) DEBUG=true
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
RNA_FILE="${DATA_DIR}/rna/CDP_merged.fa"


# switch directory according to options
cd $ROOT_DIR
if [ "$MYERS" = true ]; then
    # run tests
    for ((e=$eLow;e<=$eHigh;e+=5)); do
        for ((c=$cLow;c<=$cHigh;c+=1)); do
            OUTPUT="chr${chr}_myers_l${l}_c${c}_e${e}.out"

            echo "Running l: ${l}, c: ${c}, e: ${e}"
            echo "$EXEC_CMD -ss $RNA_FILE -ds $DNA_FILE -o $OUTPUT -i -v -l $l -fm 1 -c $c -e $e &> ${OUTPUT}.debug"
            $EXEC_CMD -ss $RNA_FILE -ds $DNA_FILE -o $OUTPUT -i -v -l $l -fm 1 -c $c -e $e &> ${OUTPUT}.debug
            #rm $OUTPUT ${OUTPUT}.summary ${OUTPUT}.debug
            rm ${OUTPUT}.debug
        done
    done
fi

