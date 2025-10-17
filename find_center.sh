#!/bin/bash

PATH_PROBE=$1
DIR_UCE=$2
OUTPUT_PATH=$3

 find $DIR_UCE -type f -name "*.fasta" | xargs -I[] bash -c \
	 "echo [] | sed 's/\.fasta//' | sed 's/.*\///' | xargs -I{} bash -c 'echo {} && grep -A3 {} $PATH_PROBE > /tmp/tmp.probe' \
	 && cat [] | sed 's/?/n/g' > /tmp/tmp.fasta \
	 && ~/anaconda3/envs/phyluce-1.7.1/bin/mafft --localpair --addfragments  /tmp/tmp.probe --ep 0.123 --keeplength --mapout /tmp/tmp.fasta \
	 1> /dev/null 2>/dev/null; tail -n \$(( ( \$(cat /tmp/tmp.probe.map |wc -l) -2) /2 )) /tmp/tmp.probe.map | head -n1" | tee ${OUTPUT_PATH}

paste <(sed -n 'p;n' $OUTPUT_PATH) <(sed -n 'n;p' $OUTPUT_PATH) > tmp
mv tmp $OUTPUT_PATH
sed -i 's/	/,/g' $OUTPUT_PATH
sed -i 's/ //g' $OUTPUT_PATH
