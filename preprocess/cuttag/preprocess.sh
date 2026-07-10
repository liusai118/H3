#!/usr/bin/env bash
set -euo pipefail

patha="/home/liusai/cut_tag/prj/syt_h3"
pathr="/home/liusai/cut_tag/prj/syt_h3/result"

path1="$patha"
path2="$pathr/1filterdata"
path3="$pathr/2bam"
path4="$pathr/3bed"
path5="$pathr/4macs2ou"
path7="$pathr/7bw"

# mkdir -p "$pathr"
mkdir -p "$path2" "$path3" "$path4" "$path5" "$path7"

while read -r id; do
  [[ -z "$id" ]] && continue

  mkdir -p "$path2/${id}" "$path3/${id}" "$path5/${id}"

  fastp \
    -i "$path1/${id}/${id}_1.fq.gz" \
    -I "$path1/${id}/${id}_2.fq.gz" \
    -o "$path2/${id}/${id}_filter_R1.fq.gz" \
    -O "$path2/${id}/${id}_filter_R2.fq.gz" \
    -h "$path2/${id}/${id}_report.html"

  # fastqc -o "$path2/${id}" -t 6 \
  #   "$path2/${id}/${id}_filter_R1.fq.gz" \
  #   "$path2/${id}/${id}_filter_R2.fq.gz"

  bowtie2 \
    --end-to-end --very-sensitive \
    --no-mixed --no-discordant \
    --phred33 -I 10 -X 700 -p 10 \
    -x "/home/liusai/index/bowtie2/mouse/ucsc/mm39" \
    -1 "$path2/${id}/${id}_filter_R1.fq.gz" \
    -2 "$path2/${id}/${id}_filter_R2.fq.gz" \
    -S "$path3/${id}/${id}_bowtie2.sam"


  samtools view -@ 10 -bS -F 0x04 \
    "$path3/${id}/${id}_bowtie2.sam" \
    | samtools sort -@ 8 \
    -o "$path3/${id}/${id}_bowtie2sort.bam"

  sambamba markdup -t 6 -r -p \
    "$path3/${id}/${id}_bowtie2sort.bam" \
    "$path3/${id}/${id}_bowtie2sort_markdup.bam"

  samtools index -@ 8 \
    "$path3/${id}/${id}_bowtie2sort_markdup.bam"


  rm "$path3/${id}/${id}_bowtie2.sam"

done < syt.txt

