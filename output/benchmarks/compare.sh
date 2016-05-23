#!/bin/bash

l=25
chr=1
eLow=10
eHigh=20
cLow=1
cHigh=3

ROOT_DIR=$PWD/..
MYERS_TMP="${ROOT_DIR}/myers.tmp"
BRUTE_TMP="${ROOT_DIR}/brute.tmp"

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

for ((e=$eLow;e<=$eHigh;e+=5)); do
    for ((c=$cLow;c<=$cHigh;c+=1)); do
        echo "Comparing l: ${l}, c: ${c}, e: ${e}..."

        OUTPUT="chr${chr}_myers_l${l}_c${c}_e${e}"
        sort ${OUTPUT}.tpx > ${MYERS_TMP} 
        OUTPUT="../mm9/original/chr${chr}_brute_l${l}_c${c}_e${e}"
        sort ${OUTPUT}.tpx > ${BRUTE_TMP}

        if [ ! -f ${MYERS_TMP} ]; then
            echo "Myers sorted output not found. Can't compare."
        elif [ ! -f ${BRUTE_TMP} ]; then
            echo "Brute sorted output not found. Can't compare."
        else
            echo "Comparing myers and brute-force... wait for failed?"
            cmp --silent ${MYERS_TMP} ${BRUTE_TMP} || echo "Test #${TEST} failed."
        fi
    
    done
done
