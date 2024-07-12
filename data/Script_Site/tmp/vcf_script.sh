#!/bin/bash

source /home/grid/miniconda3/etc/profile.d/conda.sh
conda activate genomics

mkdir -p C:\Users\aleks\OneDrive\Bureau\CHU\Test1\vcf
samtools faidx aaa
samtools index C:\Users\aleks\OneDrive\Bureau\CHU\Test1\merge_bam
bcftools mpileup -Ou -f aaa C:\Users\aleks\OneDrive\Bureau\CHU\Test1\merge_bam | bcftools call -mv -Ob -o C:\Users\aleks\OneDrive\Bureau\CHU\Test1\vcf\TTR.bcf
bcftools index C:\Users\aleks\OneDrive\Bureau\CHU\Test1\vcf\TTR.bcf
bcftools view -Oz -o C:\Users\aleks\OneDrive\Bureau\CHU\Test1\vcf\TTR.vcf.gz C:\Users\aleks\OneDrive\Bureau\CHU\Test1\vcf\TTR.bcf
tabix -p vcf C:\Users\aleks\OneDrive\Bureau\CHU\Test1\vcf\TTR.vcf.gz
gunzip -c C:\Users\aleks\OneDrive\Bureau\CHU\Test1\vcf\TTR.vcf.gz > C:\Users\aleks\OneDrive\Bureau\CHU\Test1\vcf\TTR.vcf
rm -f C:\Users\aleks\OneDrive\Bureau\CHU\Test1\vcf\TTR.bcf C:\Users\aleks\OneDrive\Bureau\CHU\Test1\vcf\TTR.vcf.gz C:\Users\aleks\OneDrive\Bureau\CHU\Test1\vcf\TTR.bcf.csi C:\Users\aleks\OneDrive\Bureau\CHU\Test1\vcf\TTR.vcf.gz.tbi
echo "Variant calling and file processing completed."
