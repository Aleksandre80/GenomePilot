#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate genomics

BASE_OUTPUT_DIR="C:\Users\aleks\OneDrive\Bureau\CHU\Test1\Basecalling"
mkdir -p "${BASE_OUTPUT_DIR}"

DORADO_BIN="/home/grid/dorado-0.7.2-linux-x64/bin/dorado"
MODEL_PATH="/home/grid/dorado-0.7.2-linux-x64/bin/dna_r10.4.1_e8.2_400bps_sup@v3.5.2"
REF_GENOME="aaa"
INPUT_DIR="C:\Users\aleks\OneDrive\Bureau\CHU\Test1\pod5"
OUTPUT_DIR="${BASE_OUTPUT_DIR}/demultiplexed_q14"
mkdir -p "${OUTPUT_DIR}"
${DORADO_BIN} basecaller -x "cuda:all" --min-qscore "14" --no-trim --emit-fastq ${MODEL_PATH} ${INPUT_DIR} | \
${DORADO_BIN} demux --kit-name "SQK-NBD114-24" --emit-fastq --output-dir "${OUTPUT_DIR}"
echo "Processing complete for C:\Users\aleks\OneDrive\Bureau\CHU\Test1\pod5 with Q-score 14"
for fastq_file in "${OUTPUT_DIR}"/*.fastq; do
    bam_file="${fastq_file%.fastq}.bam"
    echo "Aligning ${fastq_file} to reference genome..."
    minimap2 -ax map-ont "aaa" "$fastq_file" | samtools sort -o "$bam_file"
    samtools index "$bam_file"
    echo "Alignment and BAM conversion completed for ${bam_file}"
done
BASE_OUTPUT_DIR="C:\Users\aleks\OneDrive\Bureau\CHU\Test1\Basecalling"
mkdir -p "${BASE_OUTPUT_DIR}"

DORADO_BIN="/home/grid/dorado-0.7.2-linux-x64/bin/dorado"
MODEL_PATH="/home/grid/dorado-0.7.2-linux-x64/bin/dna_r10.4.1_e8.2_400bps_sup@v3.5.2"
REF_GENOME="aaa"
INPUT_DIR="C:\Users\aleks\OneDrive\Bureau\CHU\Test1\pod5"
OUTPUT_DIR="${BASE_OUTPUT_DIR}/demultiplexed_q17"
mkdir -p "${OUTPUT_DIR}"
${DORADO_BIN} basecaller -x "cuda:all" --min-qscore "17" --no-trim --emit-fastq ${MODEL_PATH} ${INPUT_DIR} | \
${DORADO_BIN} demux --kit-name "SQK-NBD114-24" --emit-fastq --output-dir "${OUTPUT_DIR}"
echo "Processing complete for C:\Users\aleks\OneDrive\Bureau\CHU\Test1\pod5 with Q-score 17"
for fastq_file in "${OUTPUT_DIR}"/*.fastq; do
    bam_file="${fastq_file%.fastq}.bam"
    echo "Aligning ${fastq_file} to reference genome..."
    minimap2 -ax map-ont "aaa" "$fastq_file" | samtools sort -o "$bam_file"
    samtools index "$bam_file"
    echo "Alignment and BAM conversion completed for ${bam_file}"
done
BASE_OUTPUT_DIR="C:\Users\aleks\OneDrive\Bureau\CHU\Test1\Basecalling"
mkdir -p "${BASE_OUTPUT_DIR}"

DORADO_BIN="/home/grid/dorado-0.7.2-linux-x64/bin/dorado"
MODEL_PATH="/home/grid/dorado-0.7.2-linux-x64/bin/dna_r10.4.1_e8.2_400bps_sup@v3.5.2"
REF_GENOME="aaa"
INPUT_DIR="C:\Users\aleks\OneDrive\Bureau\CHU\Test1\pod5"
OUTPUT_DIR="${BASE_OUTPUT_DIR}/demultiplexed_q20"
mkdir -p "${OUTPUT_DIR}"
${DORADO_BIN} basecaller -x "cuda:all" --min-qscore "20" --no-trim --emit-fastq ${MODEL_PATH} ${INPUT_DIR} | \
${DORADO_BIN} demux --kit-name "SQK-NBD114-24" --emit-fastq --output-dir "${OUTPUT_DIR}"
echo "Processing complete for C:\Users\aleks\OneDrive\Bureau\CHU\Test1\pod5 with Q-score 20"
for fastq_file in "${OUTPUT_DIR}"/*.fastq; do
    bam_file="${fastq_file%.fastq}.bam"
    echo "Aligning ${fastq_file} to reference genome..."
    minimap2 -ax map-ont "aaa" "$fastq_file" | samtools sort -o "$bam_file"
    samtools index "$bam_file"
    echo "Alignment and BAM conversion completed for ${bam_file}"
done
echo "All processes are complete."
