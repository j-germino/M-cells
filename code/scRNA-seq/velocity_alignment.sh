#!/bin/bash

echo "Fezf2 KO 1"

kb count \
    --h5ad \
    -i /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/index.idx \
    -g /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/t2g.txt \
    -x 10XV3 \
    -o Fezf2_KO_1 \
    -c1 /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/spliced_t2c.txt \
    -c2 /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/unspliced_t2c.txt \
    --workflow lamanno \
    --filter bustools \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/andersonm-Fezf_KO/Fezf_MT_S1_L001_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/andersonm-Fezf_KO/Fezf_MT_S1_L001_R2_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/andersonm-Fezf_KO/Fezf_MT_S1_L002_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/andersonm-Fezf_KO/Fezf_MT_S1_L002_R2_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/andersonm-Fezf_KO/Fezf_MT_S1_L003_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/andersonm-Fezf_KO/Fezf_MT_S1_L003_R2_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/andersonm-Fezf_KO/Fezf_MT_S1_L004_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/andersonm-Fezf_KO/Fezf_MT_S1_L004_R2_001.fastq.gz
    
echo "Fezf2 WT 1"

kb count \
    --h5ad \
    -i /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/index.idx \
    -g /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/t2g.txt \
    -x 10XV3 \
    -o Fezf2_WT_1 \
    -c1 /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/spliced_t2c.txt \
    -c2 /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/unspliced_t2c.txt \
    --workflow lamanno \
    --filter bustools \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/andersonm-Fezf_WT/Fezf_WT_S2_L001_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/andersonm-Fezf_WT/Fezf_WT_S2_L001_R2_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/andersonm-Fezf_WT/Fezf_WT_S2_L002_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/andersonm-Fezf_WT/Fezf_WT_S2_L002_R2_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/andersonm-Fezf_WT/Fezf_WT_S2_L003_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/andersonm-Fezf_WT/Fezf_WT_S2_L003_R2_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/andersonm-Fezf_WT/Fezf_WT_S2_L004_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/andersonm-Fezf_WT/Fezf_WT_S2_L004_R2_001.fastq.gz

echo "Fezf2 KO 2"

kb count \
    --h5ad \
    -i /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/index.idx \
    -g /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/t2g.txt \
    -x 10XV3 \
    -o Fezf2_KO_2 \
    -c1 /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/spliced_t2c.txt \
    -c2 /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/unspliced_t2c.txt \
    --workflow lamanno \
    --filter bustools \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/Fezf2_cKO2/andersonm-Fezf2_cKO2_S14_L001_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/Fezf2_cKO2/andersonm-Fezf2_cKO2_S14_L001_R2_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/Fezf2_cKO2/andersonm-Fezf2_cKO2_S14_L002_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/Fezf2_cKO2/andersonm-Fezf2_cKO2_S14_L002_R2_001.fastq.gz
    
echo "Fezf2 WT 2"

kb count \
    --h5ad \
    -i /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/index.idx \
    -g /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/t2g.txt \
    -x 10XV3 \
    -o Fezf2_WT_2 \
    -c1 /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/spliced_t2c.txt \
    -c2 /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/unspliced_t2c.txt \
    --workflow lamanno \
    --filter bustools \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/Fezf2_WT2/andersonm-Fezf2_WT2_S13_L001_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/Fezf2_WT2/andersonm-Fezf2_WT2_S13_L001_R2_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/Fezf2_WT2/andersonm-Fezf2_WT2_S13_L002_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Fezf2_KO/fastqs/Fezf2_WT2/andersonm-Fezf2_WT2_S13_L002_R2_001.fastq.gz
    
echo "Aire KO 1"

kb count \
    --h5ad \
    -i /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/index.idx \
    -g /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/t2g.txt \
    -x 10XV3 \
    -o Aire_KO_1 \
    -c1 /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/spliced_t2c.txt \
    -c2 /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/unspliced_t2c.txt \
    --workflow lamanno \
    --filter bustools \
    /mnt/iacchus/joe/raw_data/Thymus/Yi/Yi_RNA_fastqs/WT_KO_fastqs/Aire_KO/B6_Aire_KO_S23_L001_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Yi/Yi_RNA_fastqs/WT_KO_fastqs/Aire_KO/B6_Aire_KO_S23_L001_R2_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Yi/Yi_RNA_fastqs/WT_KO_fastqs/Aire_KO/B6_Aire_KO_S23_L002_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Yi/Yi_RNA_fastqs/WT_KO_fastqs/Aire_KO/B6_Aire_KO_S23_L002_R2_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Yi/Yi_RNA_fastqs/WT_KO_fastqs/Aire_KO/B6_Aire_KO_S23_L003_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Yi/Yi_RNA_fastqs/WT_KO_fastqs/Aire_KO/B6_Aire_KO_S23_L003_R2_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Yi/Yi_RNA_fastqs/WT_KO_fastqs/Aire_KO/B6_Aire_KO_S23_L004_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Yi/Yi_RNA_fastqs/WT_KO_fastqs/Aire_KO/B6_Aire_KO_S23_L004_R2_001.fastq.gz
    
echo "Aire WT 1"

kb count \
    --h5ad \
    -i /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/index.idx \
    -g /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/t2g.txt \
    -x 10XV3 \
    -o Aire_WT_1 \
    -c1 /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/spliced_t2c.txt \
    -c2 /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/unspliced_t2c.txt \
    --workflow lamanno \
    --filter bustools \
    /mnt/iacchus/joe/raw_data/Thymus/Yi/Yi_RNA_fastqs/WT_KO_fastqs/Aire_WT/B6_Aire_WT_S22_L001_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Yi/Yi_RNA_fastqs/WT_KO_fastqs/Aire_WT/B6_Aire_WT_S22_L001_R2_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Yi/Yi_RNA_fastqs/WT_KO_fastqs/Aire_WT/B6_Aire_WT_S22_L002_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Yi/Yi_RNA_fastqs/WT_KO_fastqs/Aire_WT/B6_Aire_WT_S22_L002_R2_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Yi/Yi_RNA_fastqs/WT_KO_fastqs/Aire_WT/B6_Aire_WT_S22_L003_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Yi/Yi_RNA_fastqs/WT_KO_fastqs/Aire_WT/B6_Aire_WT_S22_L003_R2_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Yi/Yi_RNA_fastqs/WT_KO_fastqs/Aire_WT/B6_Aire_WT_S22_L004_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Yi/Yi_RNA_fastqs/WT_KO_fastqs/Aire_WT/B6_Aire_WT_S22_L004_R2_001.fastq.gz
    
echo "Aire KO 2"

kb count \
    --h5ad \
    -i /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/index.idx \
    -g /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/t2g.txt \
    -x 10XV3 \
    -o Aire_KO_2 \
    -c1 /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/spliced_t2c.txt \
    -c2 /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/unspliced_t2c.txt \
    --workflow lamanno \
    --filter bustools \
    /mnt/iacchus/joe/raw_data/Thymus/Aire_KO/fastqs/Aire_KO-2/andersonm-02-Aire_KO2_S2_L001_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Aire_KO/fastqs/Aire_KO-2/andersonm-02-Aire_KO2_S2_L001_R2_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Aire_KO/fastqs/Aire_KO-2/andersonm-02-Aire_KO2_S2_L002_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Aire_KO/fastqs/Aire_KO-2/andersonm-02-Aire_KO2_S2_L002_R2_001.fastq.gz
    
echo "Aire WT 2"

kb count \
    --h5ad \
    -i /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/index.idx \
    -g /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/t2g.txt \
    -x 10XV3 \
    -o Aire_WT_2 \
    -c1 /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/spliced_t2c.txt \
    -c2 /mnt/iacchus/joe/reference_genomes/kallisto/mus_musculus_velocity/unspliced_t2c.txt \
    --workflow lamanno \
    --filter bustools \
    /mnt/iacchus/joe/raw_data/Thymus/Aire_KO/fastqs/Aire_WT-2/andersonm-01-Aire_WT2_S1_L001_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Aire_KO/fastqs/Aire_WT-2/andersonm-01-Aire_WT2_S1_L001_R2_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Aire_KO/fastqs/Aire_WT-2/andersonm-01-Aire_WT2_S1_L002_R1_001.fastq.gz \
    /mnt/iacchus/joe/raw_data/Thymus/Aire_KO/fastqs/Aire_WT-2/andersonm-01-Aire_WT2_S1_L002_R2_001.fastq.gz