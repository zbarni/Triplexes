#!/bin/bash

TTS_FILE='data.tts'
TFO_FILE='data.tfo'

#1
cd "output"
for TEST in {1..9}
do
    triplexator -ss "../data/tfo_${TEST}.data" -ds "../data/tts_${TEST}.data" -e 20 -c 1 -l 10 -fm 1 -i  -o "myers_${TEST}.tpx" -of 1 &> "myers_${TEST}.dbg"
    triplexator -ss "../data/tfo_${TEST}.data" -ds "../data/tts_${TEST}.data" -e 20 -c 1 -l 10 -o "brute_${TEST}.tpx" -of 1 &> "brute_${TEST}.dbg"
#    triplexator -ss "../data/tfo_${TEST}.data" -ds "../data/tts_${TEST}.data" -e 20 -c 1 -l 10 -fm 1 -i  -o "myers_${TEST}.tpx"  &> "myers_${TEST}.dbg"
#    triplexator -ss "../data/tfo_${TEST}.data" -ds "../data/tts_${TEST}.data" -e 20 -c 1 -l 10 -o "brute_${TEST}.tpx" &> "brute_${TEST}.dbg"
    sort "myers_${TEST}.tpx" > "myers.tmp"
    sort "brute_${TEST}.tpx" > "brute.tmp"
    cmp --silent "myers.tmp" "brute.tmp"  || echo "Test #${TEST} failed."
    #rm "myers.tmp" "brute.tmp"
done
