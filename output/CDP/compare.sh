#!/bin/bash
error=0
cd inverted
for f in *.out; do
    temp1="../${f}.sorted.temp1"
    temp2="../${f}.sorted.temp2"
    temp3="../${f}.sorted.temp3"
    sort -k1 -n $f > $temp1
    sort -k1 -n ../original/$f > $temp2
    sort -k1 -n ../original-local/$f > $temp3
    cmp --silent $temp1 $temp2 || echo "mismatch1 in ${f}" || error=1
    cmp --silent $temp1 $temp3 || echo "mismatch2 in ${f}" || error=1
    rm $temp1 $temp2 $temp3
done

if [ $error == 0 ]; then
    echo "No mismatches"
fi
