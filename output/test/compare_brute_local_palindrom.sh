#!/bin/bash

LARGE=""
VERBOSE=""
TEST=""
OD_LOCAL=$PWD"/results/bit_parallel_local"
OD_BIT=$PWD"/results/bit_parallel"
OD_BRUTE=$PWD"/results/brute"
OPTIONS=" -e 20 -c 2 -L -1 -l 15 -fr off "

while getopts ":vx:t:" opt; do
    case ${opt} in
        v ) 
        VERBOSE=" -of 1 "
        ;;
        x ) 
        LARGE="x"
        ;;
        t ) TEST=$OPTARG
        ;;
        \? ) echo "Usage: run.sh [-t TESTNUMBER] [-x] [-v]"
        ;;
    esac
done

if [[ $TEST != "" ]]; then
    echo "Running MYERS LOCAL #${LARGE}${TEST}"
    #valgrind --tool=memcheck triplexator -ss "data/tts_${LARGE}${TEST}.data" -ds "data/tts_${LARGE}${TEST}.data" $OPTIONS -fm 1 --bit-parallel-local  -o "${LARGE}${TEST}.tpx" -od $OD_LOCAL $VERBOSE  &> "${OD_LOCAL}/${LARGE}${TEST}.dbg"
    triplexator -ss "data/tts_${LARGE}${TEST}.data" -ds "data/tts_${LARGE}${TEST}.data" $OPTIONS -fm 1 --bit-parallel-local  -o "${LARGE}${TEST}.tpx" -od $OD_LOCAL $VERBOSE  &> "${OD_LOCAL}/${LARGE}${TEST}.dbg"

    echo "Running BRUTE FORCE #${LARGE}${TEST}"
    triplexator -ss "data/tts_${LARGE}${TEST}.data" -ds "data/tts_${LARGE}${TEST}.data" $OPTIONS -o "${LARGE}${TEST}.tpx" -od $OD_BRUTE $VERBOSE &> "${OD_BRUTE}/${LARGE}${TEST}.dbg"

    echo ""
    sort "${OD_LOCAL}/${LARGE}${TEST}.tpx" > "bpl.tmp"
    sort "${OD_BRUTE}/${LARGE}${TEST}.tpx" > "brute.tmp"
    compare_brute_local_palindrom.py $PWD/bpl.tmp $PWD/brute.tmp
#    rm bpl.tmp
#    rm brute.tmp
    exit
fi


for TEST in {1..36}
do
done
