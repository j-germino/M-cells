#!/bin/bash

computeMatrix reference-point \
    -S /mnt/iacchus/joe/raw_data/Thymus/Fezf2_chip/output/bigwigs/Fezf2_input_sorted_rmdup.bw \
       /mnt/iacchus/joe/raw_data/Thymus/Fezf2_chip/output/bigwigs/Fezf2_sorted_rmdup.bw \
    -R /home/joe/Repositories/fezf2-thymus/analysis/DE/Fezf2_DE_up_genes.bed \
       /home/joe/Repositories/fezf2-thymus/analysis/DE/Fezf2_DE_down_genes.bed \
    -p 16 \
    -b 5000 \
    -a 5000 \
    --skipZeros \
    -o /home/joe/Repositories/fezf2-thymus/analysis/ChIP_seq/de_genes_matrix.mat.gz

plotHeatmap -m /home/joe/Repositories/fezf2-thymus/analysis/ChIP_seq/de_genes_matrix.mat.gz \
    -out /home/joe/Repositories/fezf2-thymus/analysis/ChIP_seq/DE_genes_heatmap.pdf \
    --dpi 300 \
    --outFileSortedRegions /home/joe/Repositories/fezf2-thymus/analysis/ChIP_seq/DE_genes_regions_sorted.bed \
    --legendLocation upper-left \
    --colorMap Blues \

plotHeatmap -m /home/joe/Repositories/fezf2-thymus/analysis/ChIP_seq/de_genes_matrix.mat.gz \
    -out /home/joe/Repositories/fezf2-thymus/analysis/ChIP_seq/DE_genes_heatmap.png \
    --dpi 300 \
    --outFileSortedRegions /home/joe/Repositories/fezf2-thymus/analysis/ChIP_seq/DE_genes_regions_sorted.bed \
    --legendLocation upper-left \
    --colorMap Blues \