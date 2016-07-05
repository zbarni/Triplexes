#!/bin/bash

LARGE=""
VERBOSE=""
TEST=""

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
#1
cd "output"

if [[ $TEST != "" ]]; then
    echo "Running MYERS LOCAL #${LARGE}${TEST}"
    triplexator -ss "../data/tts_${LARGE}${TEST}.data" -ds "../data/tts_${LARGE}${TEST}.data" -e 20 -c 1 -L -1 -l 20 -fm 1 --bit-parallel-local  -o "local_myers_${LARGE}${TEST}.tpx" $VERBOSE  &> "local_myers_${LARGE}${TEST}.dbg"

    echo "Running BRUTE FORCE #${LARGE}${TEST}"
    triplexator -ss "../data/tts_${LARGE}${TEST}.data" -ds "../data/tts_${LARGE}${TEST}.data" -e 20 -c 1 -L -1 -l 20 -o "brute_${LARGE}${TEST}.tpx" $VERBOSE &> "brute_${LARGE}${TEST}.dbg"

    echo ""
    sort "local_myers_${LARGE}${TEST}.tpx" > "local_myers.tmp"
    sort "brute_${LARGE}${TEST}.tpx" > "brute.tmp"
    compare_brute_local_palindrom.py $PWD/local_myers.tmp $PWD/brute.tmp
    exit
fi


for TEST in {1..36}
#for TEST in {13..19}
do
    echo "Running MYERS #${LARGE}${TEST}"
    triplexator -ss "../data/tfo_${LARGE}${TEST}.data" -ds "../data/tts_${LARGE}${TEST}.data" -e 20 -c 2 -l 15 -fr off -L -1 -fm 1 --bit-parallel  -o "myers_${LARGE}${TEST}.tpx" $VERBOSE  &> "myers_${LARGE}${TEST}.dbg"

    echo "Running BRUTE FORCE #${LARGE}${TEST}"
    triplexator -ss "../data/tfo_${LARGE}${TEST}.data" -ds "../data/tts_${LARGE}${TEST}.data" -e 20 -c 2 -l 15 -fr off -L -1 -o "brute_${LARGE}${TEST}.tpx" $VERBOSE &> "brute_${LARGE}${TEST}.dbg"

    echo ""
    sort "myers_${LARGE}${TEST}.tpx" > "myers.tmp"
    sort "brute_${LARGE}${TEST}.tpx" > "brute.tmp"
    cmp --silent "myers.tmp" "brute.tmp"  || echo "==> Test #${LARGE}${TEST} failed. <=="
    #rm "myers.tmp" "brute.tmp"
done
