#!/bin/bash

TEST=""
COMP="myers"

while getopts ":t:q" opt; do
    case ${opt} in
        t ) 
        TEST=$OPTARG
        ;;
        q ) 
        COMP="qgram"
        ;;
        \? ) echo "Usage: meld [-b] [-m]"
        ;;
    esac
done

sort output/brute_${TEST}.tpx > tmp_brute
sort output/${COMP}_${TEST}.tpx > tmp_${COMP}

meld tmp_brute tmp_${COMP} &
