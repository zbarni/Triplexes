#!/bin/bash
mem=false   #memory profiling
inv="f"      #inverted if -i given, otherwise original
dna="../../../data/dna/de_genes_promoters_CDP_cDC_UP.fa"
rna="../../../data/rna/CDP_merged.fa"
l=15

while getopts ":mil:" opt; do
    case ${opt} in
        m ) mem=true
        ;;
        i ) inv="-i"
        ;;
        l ) l=$OPTARG
        ;;
        \? ) echo "Usage: run_tests.sh [-m] [-i] [-l min-length]"
        ;;
    esac
done

# switch directory according to options
if [ $inv == "-i" ]; then
    cd inverted
else
    cd original
fi

# run tests
for e in {5..15..5}; do
    for q in {1..2}; do
        for c in {1..2}; do
            output="serial_l${l}_q${q}_c${c}_e${e}.out"
            tempOutput="tempOutput"
            if [ "$mem" = true ]; then
                echo "Running memory profiler with l: ${l}, q: ${q}, c: ${c}, e: ${e}"
                valgrind --tool=massif --massif-out-file="profiling/${output}.mem" triplexator -ss $rna -ds $dna -o $tempOutput $inv -v -l $l -fm 1 -t $q -c $c -e $e 
                rm ${tempOutput}*
            else
                echo "Running l: ${l}, q: ${q}, c: ${c}, e: ${e}"
                triplexator -ss $rna -ds $dna -o $output $inv -v -l $l -fm 1 -t $q -c $c -e $e 
            fi
        done
    done
done
