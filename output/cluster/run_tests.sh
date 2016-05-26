#!/bin/bash
# IMPORTANT: must be run from output/mm9 or similar

l=25
chr=1
DATA_DIR="${PWD}/../../data"
ROOT_DIR=$PWD
EXEC_CMD="${PWD}/../../triplexator/bin/triplexator"
LOCAL=true
DEBUG=false
BRUTE=false
MYERS=false
eLow=5
eHigh=20
cLow=1
cHigh=3

while getopts ":mdbil:c:" opt; do
    case ${opt} in
        m ) MYERS=true
        ;;
        d ) DEBUG=true
        ;;
        b ) BRUTE=true
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
    cd myers
    echo "Move to directory: $PWD"
    # run tests
    for ((e=$eLow;e<=$eHigh;e+=5)); do
        for ((c=$cLow;c<=$cHigh;c+=1)); do
            OUTPUT="chr${chr}_myers_l${l}_c${c}_e${e}.LOWMEM.tpx"

            echo "Running l: ${l}, c: ${c}, e: ${e}"
            $ROOT_DIR/meta_bsub.sh bsub mm9_myers mm9_chr${chr}_l${l}_c${c}_e${e}.LOWMEM 24:00 32000 $EXEC_CMD -ss $RNA_FILE -ds $DNA_FILE -o $OUTPUT -i -v -l $l -fm 1 -c $c -e $e 

            #rm $OUTPUT ${OUTPUT}.summary ${OUTPUT}.debug
            #rm ${OUTPUT}.debug
        done
    done
fi

cd $ROOT_DIR
if [ "$BRUTE" = true ]; then
    cd brute
    echo "Move to directory: $PWD"
    # run tests
    for ((e=$eLow;e<=$eHigh;e+=5)); do
        for ((c=$cLow;c<=$cHigh;c+=1)); do
            OUTPUT="chr${chr}_brute_l${l}_c${c}_e${e}.tpx"

            echo "Running l: ${l}, c: ${c}, e: ${e}"
            $ROOT_DIR/meta_bsub.sh bsub mm9_brute mm9_chr${chr}_l${l}_c${c}_e${e} 48:00 16000 $EXEC_CMD -ss $RNA_FILE -ds $DNA_FILE -o $OUTPUT -v -l $l -c $c -e $e &

            #rm $OUTPUT ${OUTPUT}.summary
            #rm ${OUTPUT}.debug
        done
    done
fi
