#!/usr/bin/env bash

#if (( $# < 4 ))
#then
#    echo "`basename $0` [cat/bsub] name W command args"
#    exit -1
#fi
#
#command -v $4 >/dev/null 2>&1 || { echo >&2 "Command '$2' not available"; exit -1; }
#
#bin=`which $4`
args="${@:7}"

name="$2"
outputName="$3"
W="$4"
M="$5"

#bin="./src/cse_indel_corrections.py"
#args="../.local/complete.fasta ../.local/ERR142611/ERR142611.bam 4 2 -d 0 -c ch3 --only-indels -s genome_3.pickle &> /home/bz107942/costa-tmp/output/$2.dump "


$1 <<EOF
#!/usr/bin/env zsh
### Job name
#BSUB -J $name

##BSUB -P lect0008
 
### File / path where STDOUT & STDERR will be written
### %J is the job ID, %I is the array ID
#BSUB -o ${outputName}.log.%J.%I
 
### Request the time you need for execution in minutes
### The format for the parameter is: [hour:]minute,
### that means for 80 minutes you could also use this: 1:20
#BSUB -W $W
 
### Request memory you need for your job in TOTAL in MB
#BSUB -M $M
 
### Change to the work directory
C=\$LSB_JOBINDEX
cd $PWD

source /home/${USER}/.zshrc

### Execute your application
/rwthfs/rz/cluster/home/bz107942/triplexator/triplexator/bin/triplexator $args
EOF

