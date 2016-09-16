#!/bin/bash
TEST=$1
OUT=$2

cp tfo_${TEST}.data tfo_${OUT}.data
cp tts_${TEST}.data tts_${OUT}.data
vim -O tfo_${OUT}.data tts_${OUT}.data 

