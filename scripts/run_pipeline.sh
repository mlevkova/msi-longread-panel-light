#!/bin/bash

INPUT=$1
OUTPUT=$2
REF=$3

mkdir -p $OUTPUT

# Align reads
minimap2 -ax map-ont $REF $INPUT > $OUTPUT/aligned.sam

# Convert and sort
samtools view -bS $OUTPUT/aligned.sam | samtools sort -o $OUTPUT/aligned.sorted.bam

# Index
samtools index $OUTPUT/aligned.sorted.bam
