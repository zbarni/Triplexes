#!/bin/bash
if [ -f "summary.data" ]; then 
    rm "summary.data"
fi

cd inverted

for e in {5..10..1}; do
    for q in {1..5}; do
        for c in {1..3}; do
            for FILE in "chr*q${q}*c${c}*e${e}*.log"; do
                FILE=$(echo $FILE | sed 's/\.\///')
                FILE_ORIG="../original/${FILE}"
                if ([ -f $FILE ] && [ -f $FILE_ORIG ]); then
                    TRIPLEX_SEARCH=$(sed -n '70p' $FILE)
                    echo $TRIPLEX_SEARCH
                    exit
                fi

                exit 
            done
            exit
        done
    done
done

exit

for FILE in chr*.log; do
    [[ $FILE =~ q(.)_ ]] 
    Q=${BASH_REMATCH[1]}
    
    [[ $FILE =~ c(.)_ ]] 
    C=${BASH_REMATCH[1]}

    [[ $FILE =~ _e(.[0-9]*)\. ]] 
    E=${BASH_REMATCH[1]}
    echo $Q $C $E
done
