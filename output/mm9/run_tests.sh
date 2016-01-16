#!/bin/bash

#meta_bsub.sh cat mm9 mm9i1 1:00 8000 valgrind --tool=massif --massif-out-file=/home/bz107942/triplexator/output/mm9/inverted/profiling/chr1_serial_l19_q2_c1_e5.mem /home/bz107942/triplexator/triplexator/bin/triplexator -ss /home/bz107942/triplexator/data/rna/CDP_merged.fa -ds /home/bz107942/triplexator/data/dna/mm9/mm9.chr1.fa -i -v -l 19 -fm 1 -t 2 -c 1 -e 5 \> /home/bz107942/triplexator/output/mm9/inverted/chr1_serial_l19_q2_c1_e5.out

mem=false   #memory profiling
inv="f"      #inverted if -i given, otherwise original
l=19
chr=1
#dna="../../../data/dna/mm9/mm9.chr1.fa"
#rna="../../../data/rna/CDP_merged.fa"

while getopts ":mil:c:" opt; do
    case ${opt} in
        m ) mem=true
        ;;
        i ) inv="-i"
        ;;
        l ) l=$OPTARG
        ;;
        c ) chr=$OPTARG
        ;;
        \? ) echo "Usage: run_tests.sh [-m] [-i] [-l min-length]"
        ;;
    esac
done

dna="/home/bz107942/triplexator/data/dna/mm9/mm9.chr${chr}.fa"
rna="/home/bz107942/triplexator/data/rna/CDP_merged.fa"

cd /home/bz107942/triplexator/output/mm9
# switch directory according to options
if [ $inv == "-i" ]; then
    cd inverted
    echo "Move to directory: $PWD"
else
    cd original
    echo "Move to directory: $PWD"
fi

# run tests
for e in {5..5..5}; do
    for q in {4..4}; do
        for c in {1..2}; do
            output="chr${chr}_serial_l${l}_q${q}_c${c}_e${e}.out"
            tempOutput="tempOutput"

            if [ "$mem" = true ]; then
                echo "Running memory profiler with l: ${l}, q: ${q}, c: ${c}, e: ${e}"
                meta_bsub.sh bsub mm9_mem mm9_mem_chr${chr}_l${l}_q${q}_c${c}_e${e} 24:00 16000 valgrind --tool=massif --massif-out-file="/home/bz107942/triplexator/output/mm9/inverted/profiling/${output}.mem" /home/bz107942/triplexator/triplexator/bin/triplexator -ss $rna -ds $dna -o "mem_${output}" $inv -v -l $l -fm 1 -t $q -c $c -e $e 
                #rm ${tempOutput}*
            else
                echo "Running l: ${l}, q: ${q}, c: ${c}, e: ${e}"
                meta_bsub.sh bsub mm9 mm9_chr${chr}_l${l}_q${q}_c${c}_e${e} 6:00 16000 /home/bz107942/triplexator/triplexator/bin/triplexator -ss $rna -ds $dna -o $output $inv -v -l $l -fm 1 -t $q -c $c -e $e 
            fi
        done
    done
done
