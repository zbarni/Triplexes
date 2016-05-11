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
    echo "Running MYERS #${LARGE}${TEST}"
    triplexator -ss "../data/tfo_${LARGE}${TEST}.data" -ds "../data/tts_${LARGE}${TEST}.data" -e 20 -c 1 -l 25 -L -1 -fm 1 -i  -o "myers_${LARGE}${TEST}.tpx" $VERBOSE  &> "myers_${LARGE}${TEST}.dbg"

    echo "Running BRUTE FORCE #${LARGE}${TEST}"
    triplexator -ss "../data/tfo_${LARGE}${TEST}.data" -ds "../data/tts_${LARGE}${TEST}.data" -e 20 -c 1 -l 25 -L -1 -o "brute_${LARGE}${TEST}.tpx" $VERBOSE &> "brute_${LARGE}${TEST}.dbg"

    echo ""
    sort "myers_${LARGE}${TEST}.tpx" > "myers.tmp"
    sort "brute_${LARGE}${TEST}.tpx" > "brute.tmp"
    cmp --silent "myers.tmp" "brute.tmp"  || echo "==> Test #${LARGE}${TEST} failed. <=="
    exit
fi

for TEST in {3..3}
do
    echo "Running MYERS #${LARGE}${TEST}"
    triplexator -ss "../data/tfo_${LARGE}${TEST}.data" -ds "../data/tts_${LARGE}${TEST}.data" -e 20 -c 1 -l 25 -L -1 -fm 1 -i  -o "myers_${LARGE}${TEST}.tpx" $VERBOSE  &> "myers_${LARGE}${TEST}.dbg"

    echo "Running BRUTE FORCE #${LARGE}${TEST}"
    triplexator -ss "../data/tfo_${LARGE}${TEST}.data" -ds "../data/tts_${LARGE}${TEST}.data" -e 20 -c 1 -l 25 -L -1 -o "brute_${LARGE}${TEST}.tpx" $VERBOSE &> "brute_${LARGE}${TEST}.dbg"

    echo ""
    sort "myers_${LARGE}${TEST}.tpx" > "myers.tmp"
    sort "brute_${LARGE}${TEST}.tpx" > "brute.tmp"
    cmp --silent "myers.tmp" "brute.tmp"  || echo "==> Test #${LARGE}${TEST} failed. <=="
    #rm "myers.tmp" "brute.tmp"
done
