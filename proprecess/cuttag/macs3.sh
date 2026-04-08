#!/bin/bash

pathr=/home/liusai/cut_tag/prj/syt_h3/result_mm9
path_log=${pathr}/log3
path_fastq=/home/liusai/cut_tag/prj/syt_h3/result_mm9/2bam
ctrl=/home/liusai/cut_tag/prj/syt_h3/result_mm9/2bam/IgGM/IgGM_bowtie2sort.bam

mkdir -p "${pathr}" "${path_log}" "${path_fastq}"

while read -r id
do
  echo ">>> Calling narrow peaks (q=0.01) for sample: ${id}"
  macs3 callpeak \
    -t "${path_fastq}/${id}/${id}_bowtie2sort_markdup.bam" \
    -c "${ctrl}" \
    -f BAMPE \
    -g mm \
    -q 0.01 \
    --keep-dup auto \
    -n "${id}.narrow01" \
    -B \
    --outdir "${path_log}" \
    > "${path_log}/${id}.narrow01.macs3.log" 2>&1

done < syt.txt
