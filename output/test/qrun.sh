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
    echo "Running QGRAM #${LARGE}${TEST}"
    triplexator -ss "../data/tfo_${LARGE}${TEST}.data" -ds "../data/tts_${LARGE}${TEST}.data" -e 10 -c 1 -l 30 -L -1 -fm 1 -o "qgram_${LARGE}${TEST}.tpx" $VERBOSE  &> "qgram_${LARGE}${TEST}.dbg"

    echo "Running BRUTE FORCE #${LARGE}${TEST}"
    triplexator -ss "../data/tfo_${LARGE}${TEST}.data" -ds "../data/tts_${LARGE}${TEST}.data" -e 10 -c 1 -l 30 -L -1 -o "brute_${LARGE}${TEST}.tpx" $VERBOSE &> "brute_${LARGE}${TEST}.dbg"

    echo ""
    sort "qgram_${LARGE}${TEST}.tpx" > "qgram.tmp"
    sort "brute_${LARGE}${TEST}.tpx" > "brute.tmp"
    cmp --silent "qgram.tmp" "brute.tmp"  || echo "==> Test #${LARGE}${TEST} failed. <=="
    exit
fi

#for TEST in {1..11}
for TEST in {13..19}
do
    echo "Running QGRAM #${LARGE}${TEST}"
    triplexator -ss "../data/tfo_${LARGE}${TEST}.data" -ds "../data/tts_${LARGE}${TEST}.data" -e 15 -c 1 -l 30 -L -1 -fm 1   -o "qgram_${LARGE}${TEST}.tpx" $VERBOSE  &> "qgram_${LARGE}${TEST}.dbg"

    echo "Running BRUTE FORCE #${LARGE}${TEST}"
    triplexator -ss "../data/tfo_${LARGE}${TEST}.data" -ds "../data/tts_${LARGE}${TEST}.data" -e 15 -c 1 -l 30 -L -1 -o "brute_${LARGE}${TEST}.tpx" $VERBOSE &> "brute_${LARGE}${TEST}.dbg"

    echo ""
    sort "qgram_${LARGE}${TEST}.tpx" > "qgram.tmp"
    sort "brute_${LARGE}${TEST}.tpx" > "brute.tmp"
    cmp --silent "qgram.tmp" "brute.tmp"  || echo "==> Test #${LARGE}${TEST} failed. <=="
    #rm "qgram.tmp" "brute.tmp"
done
