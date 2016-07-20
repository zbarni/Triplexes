#!/bin/bash

LARGE=""
VERBOSE=""
TEST=""
OD_LOCAL=$PWD"/results/bit_parallel_local"
OD_BIT=$PWD"/results/bit_parallel"
OD_BRUTE=$PWD"/results/brute"
OPTIONS=" -e 20 -c 1 -L -1 -l 20 "

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
    #valgrind --tool=memcheck --leak-check=full    triplexator -ss "../data/tts_${LARGE}${TEST}.data" -ds "../data/tts_${LARGE}${TEST}.data" -e 20 -c 1 -L -1 -l 20 -fm 1 --bit-parallel-local  -o "local_myers_${LARGE}${TEST}.tpx" $VERBOSE  &> "local_myers_${LARGE}${TEST}.dbg"
    triplexator -ss "data/tts_${LARGE}${TEST}.data" -ds "data/tts_${LARGE}${TEST}.data" $OPTIONS -fm 1 --bit-parallel-local  -o "${LARGE}${TEST}.tpx" -od $OD_LOCAL $VERBOSE  &> "${OD_LOCAL}/${LARGE}${TEST}.dbg"

    echo "Running BRUTE FORCE #${LARGE}${TEST}"
    triplexator -ss "data/tts_${LARGE}${TEST}.data" -ds "data/tts_${LARGE}${TEST}.data" $OPTIONS -o "${LARGE}${TEST}.tpx" -od $OD_BRUTE $VERBOSE &> "${OD_BRUTE}/${LARGE}${TEST}.dbg"

    #echo "Running MYERS #${LARGE}${TEST}"
    #triplexator -ss "../data/tts_${LARGE}${TEST}.data" -ds "../data/tts_${LARGE}${TEST}.data" $OPTIONS -fm 1 --bit-parallel -o "${LARGE}${TEST}.tpx" $VERBOSE  &> "${LARGE}${TEST}.dbg"

    #echo ""
    #sort "${LARGE}${TEST}.tpx" > "myers.tmp"
    #sort "${LARGE}${TEST}.tpx" > "brute.tmp"
    #cmp --silent "myers.tmp" "brute.tmp"  || echo "==> Test #${LARGE}${TEST} failed. <=="
    ./compare_brute_local_palindrom.sh -t ${TEST}
    exit
fi


for TEST in {1..36}
do
done
