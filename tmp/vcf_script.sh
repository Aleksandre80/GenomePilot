#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate genomics

samtools faidx gh.mmi
samtools index C:\Users\aleks\OneDrive\Bureau\CHU\Test1\merge_bam
bcftools mpileup -Ou -f gh.mmi C:\Users\aleks\OneDrive\Bureau\CHU\Test1\merge_bam | bcftools call -mv -Ob -o TTR.bcf
echo "Variant calling and file processing completed."
