#!/bin/bash
cd inverted
DATA="../summary.data"
if [ -f $DATA ]; then 
    rm $DATA
fi

echo -e "qgram\tcons_error\terror_rate\tORIGINAL\tINVERTED" >> $DATA
echo "" >> $DATA

for e in {5..10..1}; do
    for q in {1..5}; do
        for c in {1..3}; do
            for FILE in $(find . -maxdepth 1 -name "chr*q${q}*c${c}*e${e}*.log"); do
                FILE=$(echo $FILE | sed 's/\.\///')
                FILE_ORIG="../original/${FILE}"
                if ([ -f $FILE ] && [ -f $FILE_ORIG ]); then
                    echo "" >> $DATA
                    [[ $( grep "triplex search only" $FILE_ORIG ) =~ ([0-9]*\.[0-9]+)\ seconds ]] 
                    TS_ORG=${BASH_REMATCH[1]}              
                    [[ $( grep "triplex search only" $FILE ) =~ ([0-9]*\.[0-9]+)\ seconds ]] 
                    TS_INV=${BASH_REMATCH[1]}              

                    [[ $( grep "IO reading" $FILE_ORIG ) =~ ([0-9]*\.[0-9]+)\ seconds ]] 
                    IO_ORG=${BASH_REMATCH[1]}              
                    [[ $( grep "IO reading" $FILE ) =~ ([0-9]*\.[0-9]+)\ seconds ]] 
                    IO_INV=${BASH_REMATCH[1]}              

                    [[ $( grep "qgram-Finder" $FILE_ORIG ) =~ ([0-9]*\.[0-9]+)\ seconds ]] 
                    QF_ORG=${BASH_REMATCH[1]}              
                    [[ $( grep "qgram-Finder" $FILE ) =~ ([0-9]*\.[0-9]+)\ seconds ]] 
                    QF_INV=${BASH_REMATCH[1]}              

                    [[ $( grep "collectSeeds" $FILE_ORIG ) =~ ([0-9]*\.[0-9]+)\ seconds ]] 
                    CS_ORG=${BASH_REMATCH[1]}              
                    [[ $( grep "collectSeeds" $FILE ) =~ ([0-9]*\.[0-9]+)\ seconds ]] 
                    CS_INV=${BASH_REMATCH[1]}              

                    echo -e "${q}\t${c}\t${e}\t${TS_ORG}\t${TS_INV}\t\ttotal search" >> $DATA
                    echo -e "${q}\t${c}\t${e}\t${IO_ORG}\t${IO_INV}\t\tio" >> $DATA
                    echo -e "${q}\t${c}\t${e}\t${CS_ORG}\t${CS_INV}\t\tcollectSeeds" >> $DATA
                    echo -e "${q}\t${c}\t${e}\t${QF_ORG}\t${QF_INV}\t\tqgramFind" >> $DATA
                fi
            done
        done
    done
done
