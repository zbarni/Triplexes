#!/bin/bash

l=19
chr=1
DATA_DIR="${PWD}/../data"
ROOT_DIR=$PWD
EXEC_CMD="./qgram"
TTS_DIR="${DATA_DIR}/tts"
TFO_DIR="${DATA_DIR}/tfo"
OUTPUT_DIR="${PWD}/../output/myers"
LABEL=""


while getopts ":mdoil:c:" opt; do
    case ${opt} in
        m ) mem=true
        ;;
        i ) INVERTED=true
        ;;
        d ) DEBUG=true
        ;;
        o ) ORIGINAL=true
        ;;
        l ) LABEL=$OPTARG
        ;;
        c ) chr=$OPTARG
        ;;
        \? ) echo "Usage: run_tests.sh [-m] [-i] [-l min-length]"
        ;;
    esac
done

# run tests
for e in {5..5..5}; do
    for q in {1..6}; do
        for c in {1..2}; do
            BASE="${OUTPUT_DIR}/chr${chr}_serial_l${l}_q${q}_c${c}_e${e}${LABEL}.out"
            TFO="chr${chr}_serial_l${l}_q${q}_c${c}_e${e}.tfo.data"
            TTS="chr${chr}_serial_l${l}_q${q}_c${c}_e${e}.tts.data"
            TFO_FILE="${TFO_DIR}/${TFO}"
            TTS_FILE="${TTS_DIR}/${TTS}"
            #echo $TFO_FILE

            #echo "${TFO_DIR}/${TFO}"
            if [[ -e $TFO_FILE && -e $TTS_FILE ]]; then 
                echo "Running l: ${l}, q: ${q}, c: ${c}, e: ${e}"
                $EXEC_CMD $TFO_FILE $TTS_FILE &> ${BASE}
            fi
        done
    done
done
