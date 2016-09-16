#!/bin/bash
error=0
cd inverted
for f in *.out; do
    temp1="../${f}.sorted.temp1"
    temp2="../${f}.sorted.temp2"
    sort -k1 -n $f > $temp1
    sort -k1 -n ../original/$f > $temp2
    cmp --silent $temp1 $temp2
    if [[ $? != 0 ]]; then
        error=1
    fi
    rm $temp1 $temp2 
done

if [[ $error == 0 ]]; then
    echo "No mismatches"
else
    echo "Mismatch. Check"
fi
