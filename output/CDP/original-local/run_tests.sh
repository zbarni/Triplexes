#!/bin/bash
l=15
for e in {5..15..5}; do
    for q in {1..2}; do
        for c in {1..2}; do
            echo "Running l: ${l}, q: ${q}, c: ${c}, e: ${e}"
            triplexator -ss ../../../data/rna/CDP_merged.fa -ds ../../../data/dna/de_genes_promoters_CDP_cDC_UP.fa -o "serial_l${l}_q${q}_c${c}_e${e}.out" -v -l $l -fm 1 -t $q -c $c -e $e #&> "serial_l${l}_q${q}_c${c}_e${e}.out"
        done
    done
done

