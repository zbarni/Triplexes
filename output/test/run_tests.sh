#!/bin/bash
# IMPORTANT: must be run from output/mm9 or similar

mem=false   #memory profiling
l=10
chr=1
DATA_DIR="${PWD}/../../data"
ROOT_DIR=$PWD
EXEC_CMD="${PWD}/../../triplexator/bin/triplexator"
LOCAL=true
DEBUG=false
ORIGINAL=false
INVERTED=false

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

DNA_FILE="${DATA_DIR}/dna/qgram2ds.fasta"
RNA_FILE="${DATA_DIR}/rna/qgram2ss.fasta"

# switch directory according to options
cd $ROOT_DIR
if [ "$INVERTED" = true ]; then
    cd inverted
    echo "Move to directory: $PWD"
fi

# run tests
for e in {5..5..5}; do
    for q in {1..1}; do
        for c in {1..2}; do
            OUTPUT="chr${chr}_serial_l${l}_q${q}_c${c}_e${e}.out"

            if [ "$mem" = true ]; then
                echo "Running memory profiler with l: ${l}, q: ${q}, c: ${c}, e: ${e}"
				if [ "$LOCAL" = true ]; then
					valgrind --tool=massif --massif-out-file="${PWD}/profiling/${OUTPUT}.mem" $EXEC_CMD -ss $RNA_FILE -ds $DNA_FILE -o "mem_${OUTPUT}" $inv -v -l $l -fm 1 -t $q -c $c -e $e 
				else
					meta_bsub.sh bsub mm9_mem mm9_mem_chr${chr}_l${l}_q${q}_c${c}_e${e} 24:00 16000 valgrind --tool=massif --massif-out-file="${PWD}/profiling/${OUTPUT}.mem" $EXEC_CMD -ss $RNA_FILE -ds $DNA_FILE -o "mem_${OUTPUT}" $inv -v -l $l -fm 1 -t $q -c $c -e $e 
				fi
            else
                echo "Running l: ${l}, q: ${q}, c: ${c}, e: ${e}"
				if [ "$LOCAL" = true ]; then
					$EXEC_CMD -ss $RNA_FILE -ds $DNA_FILE -o $OUTPUT -i -v -l $l -fm 1 -t $q -c $c -e $e &> ${OUTPUT}.debug
				else
					meta_bsub.sh bsub mm9 mm9_chr${chr}_l${l}_q${q}_c${c}_e${e} 6:00 16000 $EXEC_CMD -ss $RNA_FILE -ds $DNA_FILE -o $OUTPUT $inv -v -l $l -fm 1 -t $q -c $c -e $e 
				fi
            fi
            #rm $OUTPUT ${OUTPUT}.summary
        done
    done
done


cd $ROOT_DIR
if [ "$ORIGINAL" = true ]; then
    cd original
    echo "Move to directory: $PWD"
fi

# run tests
for e in {5..5..5}; do
    for q in {1..1}; do
        for c in {1..2}; do
            OUTPUT="chr${chr}_serial_l${l}_q${q}_c${c}_e${e}.out"

            if [ "$mem" = true ]; then
                echo "Running memory profiler with l: ${l}, q: ${q}, c: ${c}, e: ${e}"
				if [ "$LOCAL" = true ]; then
					valgrind --tool=massif --massif-out-file="${PWD}/profiling/${OUTPUT}.mem" $EXEC_CMD -ss $RNA_FILE -ds $DNA_FILE -o "mem_${OUTPUT}" $inv -v -l $l -fm 1 -t $q -c $c -e $e 
				else
					meta_bsub.sh bsub mm9_mem mm9_mem_chr${chr}_l${l}_q${q}_c${c}_e${e} 24:00 16000 valgrind --tool=massif --massif-out-file="${PWD}/profiling/${OUTPUT}.mem" $EXEC_CMD -ss $RNA_FILE -ds $DNA_FILE -o "mem_${OUTPUT}" $inv -v -l $l -fm 1 -t $q -c $c -e $e 
				fi
            else
                echo "Running l: ${l}, q: ${q}, c: ${c}, e: ${e}"
				if [ "$LOCAL" = true ]; then
					$EXEC_CMD -ss $RNA_FILE -ds $DNA_FILE -o $OUTPUT -v -l $l -fm 1 -t $q -c $c -e $e &> ${OUTPUT}.debug
				else
					meta_bsub.sh bsub mm9 mm9_chr${chr}_l${l}_q${q}_c${c}_e${e} 6:00 16000 $EXEC_CMD -ss $RNA_FILE -ds $DNA_FILE -o $OUTPUT $inv -v -l $l -fm 1 -t $q -c $c -e $e 
				fi
            fi
            #rm $OUTPUT ${OUTPUT}.summary
        done
    done
done
