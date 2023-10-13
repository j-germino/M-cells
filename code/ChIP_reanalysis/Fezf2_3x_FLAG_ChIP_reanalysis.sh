cd /mnt/iacchus/joe/raw_data/Thymus/Fezf2_ChIP_seq/fastqs
mkdir ../bowtie_output

echo Fezf2_flag_rep_1
# Align and save sorted BAM file
bowtie2 -p 16 -q --local -U \
    Fezf2_flag_rep_1.fastq.gz -x /mnt/iacchus/joe/reference_genomes/bowtie2/mm10/mm10 \
    | samtools view -hb -t 16 \
    | samtools sort -t 16 > ../bowtie_output/Fezf2_flag_rep_1.bam
# Filter on MAPQ > 20
samtools view -q 20 -b -t 16 \
    -o ../bowtie_output/Fezf2_flag_rep_1_filtered.bam \
    ../bowtie_output/Fezf2_flag_rep_1.bam
# Mark duplicate with picard
java -jar ~/bin/picard.jar MarkDuplicates \
    I=../bowtie_output/Fezf2_flag_rep_1_filtered.bam \
    O=../bowtie_output/Fezf2_flag_rep_1_filtered_marked_dups.bam \
    M=../bowtie_output/Fezf2_flag_rep_1_filtered_marked_dup_metrics.txt \
    ASSUME_SORTED=true \
    REMOVE_DUPLICATES=true
# Cleanup by removing intermediate BAM file
rm ../bowtie_output/Fezf2_flag_rep_1_filtered.bam
# Filter blacklist regions
bedtools intersect \
    -abam ../bowtie_output/Fezf2_flag_rep_1_filtered_marked_dups.bam \
    -b /mnt/iacchus/joe/reference_genomes/blacklists/mm10-blacklist.v2.bed \
    -v > ../bowtie_output/tmp.bam
mv ../bowtie_output/tmp.bam ../bowtie_output/Fezf2_flag_rep_1_filtered_marked_dups.bam
# Create the bam index
samtools index ../bowtie_output/Fezf2_flag_rep_1_filtered_marked_dups.bam

echo Fezf2_flag_rep_2
# Align and save sorted BAM file
bowtie2 -p 16 -q --local \
    -U Fezf2_flag_rep_2.fastq.gz -x /mnt/iacchus/joe/reference_genomes/bowtie2/mm10/mm10 \
    | samtools view -hb -t 16 \
    | samtools sort -t 16 > ../bowtie_output/Fezf2_flag_rep_2.bam
# Filter on MAPQ > 20
samtools view -q 20 -b -t 16 \
    -o ../bowtie_output/Fezf2_flag_rep_2_filtered.bam ../bowtie_output/Fezf2_flag_rep_2.bam
# Mark duplicate with picard
java -jar ~/bin/picard.jar MarkDuplicates \
    I=../bowtie_output/Fezf2_flag_rep_2_filtered.bam \
    O=../bowtie_output/Fezf2_flag_rep_2_filtered_marked_dups.bam \
    M=../bowtie_output/Fezf2_flag_rep_2_filtered_marked_dup_metrics.txt \
    ASSUME_SORTED=true \
    REMOVE_DUPLICATES=true
# Cleanup by removing intermediate BAM file
rm ../bowtie_output/Fezf2_flag_rep_2_filtered.bam
# Filter blacklist regions
bedtools intersect \
    -abam ../bowtie_output/Fezf2_flag_rep_2_filtered_marked_dups.bam \
    -b /mnt/iacchus/joe/reference_genomes/blacklists/mm10-blacklist.v2.bed \
    -v > ../bowtie_output/tmp.bam
mv ../bowtie_output/tmp.bam ../bowtie_output/Fezf2_flag_rep_2_filtered_marked_dups.bam
# Create the bam index
samtools index ../bowtie_output/Fezf2_flag_rep_2_filtered_marked_dups.bam 

echo IgG_rep_1
# Align and save sorted BAM file
bowtie2 -p 16 -q --local \
    -U IgG_rep_1.fastq.gz -x /mnt/iacchus/joe/reference_genomes/bowtie2/mm10/mm10 \
    | samtools view -hb -t 16 \
    | samtools sort -t 16 > ../bowtie_output/IgG_rep_1.bam
# Filter on MAPQ > 20
samtools view -q 20 -b -t 16 \
    -o ../bowtie_output/IgG_rep_1_filtered.bam ../bowtie_output/IgG_rep_1.bam
# Mark duplicate with picard
java -jar ~/bin/picard.jar MarkDuplicates \
    I=../bowtie_output/IgG_rep_1_filtered.bam \
    O=../bowtie_output/IgG_rep_1_filtered_marked_dups.bam \
    M=../bowtie_output/IgG_rep_1_filtered_marked_dup_metrics.txt \
    ASSUME_SORTED=true \
    REMOVE_DUPLICATES=true
# Cleanup by removing intermediate BAM file
rm ../bowtie_output/IgG_rep_1_filtered.bam
# Filter blacklist regions
bedtools intersect \
    -abam ../bowtie_output/IgG_rep_1_filtered_marked_dups.bam \
    -b /mnt/iacchus/joe/reference_genomes/blacklists/mm10-blacklist.v2.bed \
    -v > ../bowtie_output/tmp.bam
mv ../bowtie_output/tmp.bam ../bowtie_output/IgG_rep_1_filtered_marked_dups.bam
# Create the bam index
samtools index ../bowtie_output/IgG_rep_1_filtered_marked_dups.bam

echo IgG_rep_2
# Align and save sorted BAM file
bowtie2 -p 16 -q --local \
    -U IgG_rep_2.fastq.gz -x /mnt/iacchus/joe/reference_genomes/bowtie2/mm10/mm10 \
    | samtools view -hb -t 16 \
    | samtools sort -t 16 > ../bowtie_output/IgG_rep_2.bam
# Filter on MAPQ > 20
samtools view -q 20 -b -t 16 \
    -o ../bowtie_output/IgG_rep_2_filtered.bam ../bowtie_output/IgG_rep_2.bam
# Mark duplicate with picard
java -jar ~/bin/picard.jar MarkDuplicates \
    I=../bowtie_output/IgG_rep_2_filtered.bam \
    O=../bowtie_output/IgG_rep_2_filtered_marked_dups.bam \
    M=../bowtie_output/IgG_rep_2_filtered_marked_dup_metrics.txt \
    ASSUME_SORTED=true \
    REMOVE_DUPLICATES=true
# Cleanup by removing intermediate BAM file
rm ../bowtie_output/IgG_rep_2_filtered.bam
# Filter blacklist regions
bedtools intersect \
    -abam ../bowtie_output/IgG_rep_2_filtered_marked_dups.bam \
    -b /mnt/iacchus/joe/reference_genomes/blacklists/mm10-blacklist.v2.bed \
    -v > ../bowtie_output/tmp.bam
mv ../bowtie_output/tmp.bam ../bowtie_output/IgG_rep_2_filtered_marked_dups.bam
# Create the bam index
samtools index ../bowtie_output/IgG_rep_2_filtered_marked_dups.bam

echo input_rep_1
# Align and save sorted BAM file
bowtie2 -p 16 -q --local -U input_rep_1.fastq.gz -x /mnt/iacchus/joe/reference_genomes/bowtie2/mm10/mm10 \
    | samtools view -hb -t 16 \
    | samtools sort -t 16 > ../bowtie_output/input_rep_1.bam
# Filter on MAPQ > 20
samtools view -q 20 -b -t 16 \
    -o ../bowtie_output/input_rep_1_filtered.bam ../bowtie_output/input_rep_1.bam
# Mark duplicate with picard
java -jar ~/bin/picard.jar MarkDuplicates \
    I=../bowtie_output/input_rep_1_filtered.bam \
    O=../bowtie_output/input_rep_1_filtered_marked_dups.bam \
    M=../bowtie_output/input_rep_1_filtered_marked_dup_metrics.txt \
    ASSUME_SORTED=true \
    REMOVE_DUPLICATES=true
# Cleanup by removing intermediate BAM file
rm ../bowtie_output/input_rep_1_filtered.bam
# Filter blacklist regions
bedtools intersect \
    -abam ../bowtie_output/input_rep_1_filtered_marked_dups.bam \
    -b /mnt/iacchus/joe/reference_genomes/blacklists/mm10-blacklist.v2.bed \
    -v > ../bowtie_output/tmp.bam
mv ../bowtie_output/tmp.bam ../bowtie_output/input_rep_1_filtered_marked_dups.bam
# Create the bam index
samtools index ../bowtie_output/input_rep_1_filtered_marked_dups.bam

echo input_rep_2
# Align and save sorted BAM file
bowtie2 -p 16 -q --local \
    -U input_rep_2.fastq.gz -x /mnt/iacchus/joe/reference_genomes/bowtie2/mm10/mm10 \
    | samtools view -hb -t 16 \
    | samtools sort -t 16 > ../bowtie_output/input_rep_2.bam
# Filter on MAPQ > 20
samtools view -q 20 -b -t 16 \
    -o ../bowtie_output/input_rep_2_filtered.bam ../bowtie_output/input_rep_2.bam
# Mark duplicate with picard
java -jar ~/bin/picard.jar MarkDuplicates \
    I=../bowtie_output/input_rep_2_filtered.bam \
    O=../bowtie_output/input_rep_2_filtered_marked_dups.bam \
    M=../bowtie_output/input_rep_2_filtered_marked_dup_metrics.txt \
    ASSUME_SORTED=true \
    REMOVE_DUPLICATES=true
# Cleanup by removing intermediate BAM file
rm ../bowtie_output/input_rep_2_filtered.bam
# Filter blacklist regions
bedtools intersect \
    -abam ../bowtie_output/input_rep_2_filtered_marked_dups.bam \
    -b /mnt/iacchus/joe/reference_genomes/blacklists/mm10-blacklist.v2.bed \
    -v > ../bowtie_output/tmp.bam
mv ../bowtie_output/tmp.bam ../bowtie_output/input_rep_2_filtered_marked_dups.bam
# Create the bam index
samtools index ../bowtie_output/input_rep_2_filtered_marked_dups.bam

# Peak calling (correct way)
mkdir ../MACS2_output
# Fezf2-FLAG to input
mkdir ../MACS2_output/Fezf2_input
macs2 callpeak \
    -t ../bowtie_output/Fezf2_flag_rep_1_filtered_marked_dups.bam ../bowtie_output/Fezf2_flag_rep_2_filtered_marked_dups.bam \
	-c ../bowtie_output/input_rep_1_filtered_marked_dups.bam ../bowtie_output/input_rep_2_filtered_marked_dups.bam \
 	-f BAM \
    -g mm \
	-n Fezf2-input \
    --keep-dup all \
    --call-summits \
	--outdir ../MACS2_output/Fezf2_input

# Fezf2-FLAG to IgG
mkdir ../MACS2_output/Fezf2_IgG
macs2 callpeak \
    -t ../bowtie_output/Fezf2_flag_rep_1_filtered_marked_dups.bam ../bowtie_output/Fezf2_flag_rep_2_filtered_marked_dups.bam \
	-c ../bowtie_output/IgG_rep_1_filtered_marked_dups.bam ../bowtie_output/IgG_rep_2_filtered_marked_dups.bam \
 	-f BAM \
    -g mm \
	-n Fezf2-IgG \
    --keep-dup all \
    --call-summits \
	--outdir ../MACS2_output/Fezf2_IgG

# Peak calling (published parameters)
# Fezf2-FLAG to input
mkdir ../MACS2_output/Fezf2_input_repeat
macs2 callpeak \
    -t ../bowtie_output/Fezf2_flag_rep_1_filtered_marked_dups.bam ../bowtie_output/Fezf2_flag_rep_2_filtered_marked_dups.bam \
	-c ../bowtie_output/input_rep_1_filtered_marked_dups.bam ../bowtie_output/input_rep_2_filtered_marked_dups.bam \
 	-f BAM \
    -g mm \
	-n Fezf2-input-repeat \
    --nomodel \
    --nolambda \
    --keep-dup all \
    --call-summits \
	--outdir ../MACS2_output/Fezf2_input_repeat

# Fezf2-FLAG to IgG
mkdir ../MACS2_output/Fezf2_IgG_repeat
macs2 callpeak \
    -t ../bowtie_output/Fezf2_flag_rep_1_filtered_marked_dups.bam ../bowtie_output/Fezf2_flag_rep_2_filtered_marked_dups.bam \
	-c ../bowtie_output/IgG_rep_1_filtered_marked_dups.bam ../bowtie_output/IgG_rep_2_filtered_marked_dups.bam \
 	-f BAM \
    -g mm \
	-n Fezf2-IgG-repeat \
    --nomodel \
    --nolambda \
    --keep-dup all \
    --call-summits \
	--outdir ../MACS2_output/Fezf2_IgG_repeat