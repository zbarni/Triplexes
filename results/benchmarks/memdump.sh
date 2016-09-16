#!/bin/bash

while true; do
    PID=`pgrep "triplexator"`
    if [[ $PID != "" ]] #pgrep "triplexator" > /dev/null
    then
        LINE=`sudo ps_mem -p $PID | head -n 5 | tail -1`
        echo -e `date` ' \t '  $LINE
        sleep 3
    else
        break;
    fi
done
