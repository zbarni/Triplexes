#!/bin/bash
# IMPORTANT: must be run from output/mm9 or similar

mem=false   #memory profiling
l=19
chr=1
DATA_DIR="${PWD}/../../data"
ROOT_DIR=$PWD
EXEC_CMD="${PWD}/../../triplexator/bin/triplexator"
LOCAL=true
DEBUG=false
ORIGINAL=false
INVERTED=false
qLow=2
qHigh=7
eLow=5
eHigh=10
cLow=1
cHigh=2

while getopts ":mdoil:c:" opt; do
    case ${opt} in
        m ) mem=true
        ;;
        i ) INVERTED=true
        ;;
        d ) DEBUG=true
        ;;
        o ) ORIGINAL=true
        ;;
        l ) l=$OPTARG
        ;;
        c ) chr=$OPTARG
        ;;
        \? ) echo "Usage: run_tests.sh [-m] [-i] [-l min-length]"
        ;;
    esac
done

DNA_FILE="${DATA_DIR}/dna/mm9/mm9.chr${chr}.fa"
RNA_FILE="${DATA_DIR}/rna/CDP_merged.fa"


# switch directory according to options
cd $ROOT_DIR
if [ "$INVERTED" = true ]; then
    cd data_inverted
    echo "Move to directory: $PWD"
    # run tests
    for ((e=$eLow;e<=$eHigh;e+=5)); do
        for ((q=$qLow;q<=$qHigh;q+=2)); do
            for ((c=$cLow;c<=$cHigh;c+=1)); do
                OUTPUT="chr${chr}_serial_l${l}_q${q}_c${c}_e${e}.out"

                if [ "$mem" = true ]; then
                    echo "Running callgrind profiler with l: ${l}, q: ${q}, c: ${c}, e: ${e}"
                    if [ "$LOCAL" = true ]; then
                        valgrind --tool=callgrind --callgrind-out-file="${PWD}/profiling/${OUTPUT}.call" $EXEC_CMD -ss $RNA_FILE -ds $DNA_FILE -o "call_${OUTPUT}" -i -v -l $l -fm 1 -t $q -c $c -e $e 
                    else
                        meta_bsub.sh bsub mm9_mem mm9_mem_chr${chr}_l${l}_q${q}_c${c}_e${e} 24:00 16000 valgrind --tool=callgrind --callgrind-out-file="${PWD}/profiling/${OUTPUT}.mem" $EXEC_CMD -ss $RNA_FILE -ds $DNA_FILE -o "mem_${OUTPUT}" $inv -v -l $l -fm 1 -t $q -c $c -e $e 
                    fi
                else
                    echo "Running l: ${l}, q: ${q}, c: ${c}, e: ${e}"
                    if [ "$LOCAL" = true ]; then
                        $EXEC_CMD -ss $RNA_FILE -ds $DNA_FILE -o $OUTPUT -i -v -l $l -fm 1 -t $q -c $c -e $e &> ${OUTPUT}.data
                    else
                        meta_bsub.sh bsub mm9 mm9_chr${chr}_l${l}_q${q}_c${c}_e${e} 6:00 16000 $EXEC_CMD -ss $RNA_FILE -ds $DNA_FILE -o $OUTPUT -i -v -l $l -fm 1 -t $q -c $c -e $e 
                    fi
                fi
                rm $OUTPUT ${OUTPUT}.summary ${OUTPUT}.log
            done
        done
    done
fi


cd $ROOT_DIR
if [ "$ORIGINAL" = true ]; then
    cd data_original
    echo "Move to directory: $PWD"
# run tests
    for ((e=$eLow;e<=$eHigh;e+=5)); do
        for ((q=$qLow;q<=$qHigh;q+=2)); do
            for ((c=$cLow;c<=$cHigh;c+=1)); do
                OUTPUT="chr${chr}_serial_l${l}_q${q}_c${c}_e${e}.out"

                if [ "$mem" = true ]; then
                    echo "Running callgrind profiler with l: ${l}, q: ${q}, c: ${c}, e: ${e}"
                    if [ "$LOCAL" = true ]; then
                        valgrind --tool=callgrind --callgrind-out-file="${PWD}/profiling/${OUTPUT}.call" $EXEC_CMD -ss $RNA_FILE -ds $DNA_FILE -o "call_${OUTPUT}" -v -l $l -fm 1 -t $q -c $c -e $e 
                    else
                        meta_bsub.sh bsub mm9_mem mm9_mem_chr${chr}_l${l}_q${q}_c${c}_e${e} 24:00 16000 valgrind --tool=callgrind --callgrind-out-file="${PWD}/profiling/${OUTPUT}.mem" $EXEC_CMD -ss $RNA_FILE -ds $DNA_FILE -o "mem_${OUTPUT}" $inv -v -l $l -fm 1 -t $q -c $c -e $e 
                    fi
                else
                    echo "Running l: ${l}, q: ${q}, c: ${c}, e: ${e}"
                    if [ "$LOCAL" = true ]; then
                        $EXEC_CMD -ss $RNA_FILE -ds $DNA_FILE -o $OUTPUT -v -l $l -fm 1 -t $q -c $c -e $e &> ${OUTPUT}.data
                    else
                        meta_bsub.sh bsub mm9 mm9_chr${chr}_l${l}_q${q}_c${c}_e${e} 6:00 16000 $EXEC_CMD -ss $RNA_FILE -ds $DNA_FILE -o $OUTPUT -v -l $l -fm 1 -t $q -c $c -e $e 
                    fi
                fi
                rm $OUTPUT ${OUTPUT}.summary ${OUTPUT}.log
            done
        done
    done
fi
