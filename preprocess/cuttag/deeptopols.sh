#!/usr/bin/env bash
set -euo pipefail

patha="/home/liusai/cut_tag/prj/syt_h3"
pathr="/home/liusai/cut_tag/prj/syt_h3/result"

path3="$pathr/2bam"
path5="$pathr/4macs2ou"
path6="$pathr/6matrix"
path7="$pathr/7bw"

mkdir -p "$path3" "$path5" "$path6" "$path7"
tss_bed="/home/liusai/index/reference/mouse/ucsc/mm9/mm9.bed"
if [[ ! -f "$tss_bed" ]]; then
  echo "ERROR: Missing TSS bed: $tss_bed"
  exit 1
fi
while IFS= read -r id || [[ -n "$id" ]]; do
  id="${id//$'\r'/}"
  [[ -z "$id" ]] && continue
  [[ "$id" =~ ^# ]] && continue

  bam1="${path3}/${id}/${id}_bowtie2sort_markdup.bam"
  if [[ ! -f "$bam1" ]]; then
    echo "Missing $bam1"
    continue
  fi

  bw="${path7}/${id}.bw"

  bamCoverage -b "$bam1" \
    -of bigwig \
    -o "$bw" \
    -p 20 \
    --ignoreDuplicates \
    --binSize 10 \
    --normalizeUsing RPKM

  matrix="${path6}/${id}_matrix_TSS.gz"
  regions="${path6}/${id}_regions_test_genes.bed"

  computeMatrix reference-point \
    --missingDataAsZero \
    --referencePoint TSS \
    -p 15 \
    -b 3000 -a 3000 \
    -R "$tss_bed" \
    -S "$bw" \
    --skipZeros \
    -o "$matrix" \
    --outFileSortedRegions "$regions"

  plotHeatmap -m "$matrix" \
    -out "${path6}/${id}_Heatmap.png" \
    --colorMap RdYlBu_r

  plotHeatmap -m "$matrix" \
    -out "${path6}/${id}_Heatmap.pdf" \
    --plotFileFormat pdf \
    --dpi 720 \
    --colorMap RdYlBu_r

  plotProfile -m "$matrix" \
    -out "${path6}/${id}_Profile.png"

  plotProfile -m "$matrix" \
    -out "${path6}/${id}_Profile.pdf" \
    --plotFileFormat pdf \
    --perGroup \
    --dpi 720

done < "syt.txt"
