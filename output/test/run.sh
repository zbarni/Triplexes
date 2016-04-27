#!/bin/bash

TTS_FILE='data.tts'
TFO_FILE='data.tfo'

#1
cd "output"
for TEST in {1..4}
do
    triplexator -ss "../data/tfo_${TEST}.data" -ds "../data/tts_${TEST}.data" -e 10 -c 1 -l 10 -fm 1 -i  -o "myers_${TEST}.tpx" -of 1 &> "myers_${TEST}.dbg"
    triplexator -ss "../data/tfo_${TEST}.data" -ds "../data/tts_${TEST}.data" -e 10 -c 1 -l 10 -o "brute_${TEST}.tpx" -of 1 &> "brute_${TEST}.dbg"
done
