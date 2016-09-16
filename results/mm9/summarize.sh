#!/bin/bash
# IMPORTANT: must be run from output/mm9 or similar

mem=false   #memory profiling
inv="f"      #inverted if -i given, otherwise original
l=19
chr=1
DATA_DIR="${PWD}/../../data"
ROOT_DIR=$PWD
EXEC_CMD="${PWD}/../../triplexator/bin/triplexator"
LOCAL=true

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

DNA_FILE="${DATA_DIR}/dna/mm9/mm9.chr${chr}.fa"
RNA_FILE="${DATA_DIR}/rna/CDP_merged.fa"

# switch directory according to options
if [ $inv == "-i" ]; then
    cd inverted/chr$chr
    echo "Move to directory: $PWD"
else
    cd original/chr$chr
    echo "Move to directory: $PWD"
fi

IFS=' '
for f in ./*.log; do
	TOTAL=$(grep -R "Finished program within" $f | cut -d " " -f 7)
	IO=$(grep -R "Time for ds IO reading" $f | cut -d " " -f 10)
	SEARCH=$(grep -R "triplex search only" $f | cut -d " " -f 9)
	IOP=$(bc -l <<< "($IO / $TOTAL) * 100")
	SEARCHP=$(bc -l <<< "($SEARCH / $TOTAL) * 100")
	echo $f >> summary.data
	echo -e "\ttotal\t\t"$TOTAL >> summary.data
	echo -e "\tio\t\t\t"$IO >> summary.data
	echo -e "\tio %\t\t"$IOP >> summary.data
	echo -e "\tsearch\t\t"$SEARCH >> summary.data
	echo -e "\tsearch %\t"$SEARCHP >> summary.data
	echo "" >> summary.data
done
