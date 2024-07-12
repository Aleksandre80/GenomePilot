#!/bin/bash

mkdir -p "C:\Users\aleks\OneDrive\Bureau\CHU\Test1\merge_bam"
samtools merge "C:\Users\aleks\OneDrive\Bureau\CHU\Test1\merge_bam/merged.bam" "C:\Users\aleks\OneDrive\Bureau\CHU\Test1\pod5"/*.bam
echo "Merging complete for BAM files in C:\Users\aleks\OneDrive\Bureau\CHU\Test1\pod5"

