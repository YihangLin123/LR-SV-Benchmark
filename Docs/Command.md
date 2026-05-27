# Complete SV Analysis Pipeline Command Manual

## Table of Contents

1. [SV Simulate](https://www.doubao.cn)

2. [Read Alignment](https://www.doubao.cn)

3. [Downsample Alignment File](https://www.doubao.cn)

4. [SV Call](https://www.doubao.cn)

5. [Somatic SV Filtering Strategies](https://www.doubao.cn)

6. [SV Evaluation](https://www.doubao.cn)

## 1\. SV Simulate

### Tool Versions

VISOR \(v1\.1\.2\), hapdiff \(v0\.9\)

### 1\.1 Germline SV Simulation

\# Simulate five types of germline structural variations \(DEL/INS/INV/DUP/TRA\)

```bash
# Deletion simulation
VISOR HACk -g reference_genome.fa -b sim_del.bed -o donor_genome_del
# Insertion simulation
VISOR HACk -g reference_genome.fa -b sim_ins.bed -o donor_genome_ins
# Inversion simulation
VISOR HACk -g reference_genome.fa -b sim_inv.bed -o donor_genome_inv
# Duplication simulation
VISOR HACk -g reference_genome.fa -b sim_dup.bed -o donor_genome_dup
# Translocation simulation
VISOR HACk -g reference_genome.fa -b sim_tra.bed -o donor_genome_tra
```

\# Generate 30X sequencing data for different long\-read platforms

```bash
# PacBio CLR simulation
VISOR LASeR -g reference_genome.fa -s donor_genome_del -b LASeR.bed -o data_clr_del_30x --coverage 30 --threads 32 --read_type pacbio --error_model pacbio2016 --qscore_model pacbio2016 --tag

# PacBio CCS simulation
VISOR LASeR -g reference_genome.fa -s donor_genome_del -b LASeR.bed -o data_ccs_del_30x --coverage 30 --threads 32 --read_type pacbio --error_model pacbio2021 --qscore_model pacbio2021 --tag

# ONT simulation
VISOR LASeR -g reference_genome.fa -s donor_genome_del -b LASeR.bed -o data_ont_del_30x --coverage 30 --threads 32 --read_type nanopore --error_model nanopore2023 --qscore_model nanopore2023 --tag
```

### 1\.2 Somatic SV Simulation

```bash
# Generate somatic variation background using paternal and maternal haplotypes
python hapdiff.py --reference reference_genome.fa --pat test1.fasta --mat test2.fasta --out-dir output --threads 32

# Merge bam files to construct tumor sequencing alignment file
samtools merge -@ 32 CHM1_CHM13_Tumor.bam CHM1.bam CHM13.bam
```

---

## 2\. Read Alignment

### Tool Versions

minimap2 \(v2\.30\), pbmm2 \(v1\.17\.0\), NGMLR \(v0\.2\.7\), samtools \(v1\.22\.1\)

### 2\.1 minimap2 Alignment

```bash
# PacBio CLR alignment
minimap2 -ax map-pb -t 32 --MD -R '@RG\tID:sampleId\tSM:sampleId' reference_genome.fa test.fastq > test.sam
# PacBio CCS/HIFI alignment
minimap2 -ax map-hifi -t 32 --MD -R '@RG\tID:sampleId\tSM:sampleId' reference_genome.fa test.fastq > test.sam
# ONT alignment
minimap2 -ax map-ont -t 32 --MD -R '@RG\tID:sampleId\tSM:sampleId' reference_genome.fa test.fastq > test.sam

# SAM to BAM conversion, sorting and indexing
samtools view -@ 32 -u test.sam | samtools sort -@ 32 -T test.bam.sort-tmp > test.bam
samtools index test.bam
```

### 2\.2 pbmm2 Alignment \(PacBio Specific\)

```bash
# CLR preset alignment
pbmm2 align reference_genome.fa test.fastq test.bam --preset SUBREAD --rg '@RG\tID:sampleId\tSM:sampleId'
# CCS preset alignment
pbmm2 align reference_genome.fa test.fastq test.bam --preset CCS --rg '@RG\tID:sampleId\tSM:sampleId'
# ONT alignment with CCS preset
pbmm2 align reference_genome.fa test.fastq test.bam --preset CCS --rg '@RG\tID:sampleId\tSM:sampleId'

# Standard sorting and indexing
samtools view -@ 32 -u test.sam | samtools sort -@ 32 -T test.bam.sort-tmp > test.bam
samtools index test.bam
```

### 2\.3 NGMLR Long\-read Alignment

```bash
# PacBio CLR alignment
ngmlr -t 32 --rg-id sampleId --rg-sm sampleId --rg-pl Pacbio -r reference_genome.fa -q test.fastq -o test.sam -x pacbio
# PacBio CCS alignment
ngmlr -t 32 --rg-id sampleId --rg-sm sampleId --rg-pl Pacbio -r reference_genome.fa -q test.fastq -o test.sam -x pacbio
# ONT alignment
ngmlr -t 32 --rg-id sampleId --rg-sm sampleId --rg-pl Nanopore -r reference_genome.fa -q test.fastq -o test.sam -x ont

# Sort and index aligned bam file
samtools view -@ 32 -u test.sam | samtools sort -@ 32 -T test.bam.sort-tmp > test.bam
samtools index test.bam
```

---

## 3\. Downsample Alignment File

### Tool Version

samtools \(v1\.22\.1\)

```bash
# Downsample bam file to 66% coverage
samtools view -h -s 0.66 test.bam | samtools view -bS > Downsample_test.bam
```

---

## 4\. SV Call

### Tool List \&amp; Versions

cuteSV2 \(v2\.1\.3\), Sniffles2 \(v2\.6\.3\), Debreak \(v1\.0\.2\), SVIM \(v2\.0\.0\), SVDF, SVHunter \(v1\.1\.0\), NanoVar \(v1\.8\.3\), SVision \(v1\.4\.0\), SVision\-pro \(v2\.5\.0\), Dysgu \(v1\.7\.0\), Delly2 \(v1\.5\.0\), SAVANA \(v1\.6\.0\), nanomonsv \(v0\.8\.0\), Severus \(v1\.6\.0\)

### 4\.1 cuteSV2

#### Germline SV Calling

```bash
# CLR data
cuteSV test.bam reference_genome.fa output/test.vcf work_tmp --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5 -t 16 -l 30 -s min_supporting_reads --genotype
# CCS data
cuteSV test.bam reference_genome.fa output/test.vcf work_tmp --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 -t 16 -l 30 -s min_supporting_reads --genotype
# ONT data
cuteSV test.bam reference_genome.fa output/test.vcf work_tmp --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 -t 16 -l 30 -s min_supporting_reads --genotype
```

#### Somatic SV Calling

```bash
# CCS tumor data
cuteSV tumor.bam reference_genome.fa output/tumor.vcf work_tmp --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 -t 16 -l 50 -s min_supporting_reads
# ONT tumor data
cuteSV tumor.bam reference_genome.fa output/tumor.vcf work_tmp --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 -t 16 -l 50 -s min_supporting_reads
```

### 4\.2 Sniffles2

#### Germline SV Calling

```bash
sniffles -i test.bam -v output/test.vcf -t 16 --minsvlen 30
```

#### Somatic SV Calling

```bash
# Tumor sample
sniffles --input tumor.bam --vcf output/tumor.vcf --snf output/tumor.snf --tandem-repeats TANDEM_REPEATS --reference reference_genome.fa -t 16 --minsvlen 50 --allow-overwrite
# Normal control sample
sniffles --input normal.bam --vcf output/normal.vcf --snf output/normal.snf --tandem-repeats TANDEM_REPEATS --reference reference_genome.fa -t 16 --minsvlen 50 --allow-overwrite
```

### 4\.3 Debreak

#### Germline SV Calling

```bash
debreak --bam test.bam -o output -m min_supporting_reads --min_size 30 --ref reference_genome.fa -t 16
```

#### Somatic SV Calling

```bash
# Tumor sample
debreak --bam tumor.bam -o output_tumor -m min_supporting_reads --min_size 50 --ref reference_genome.fa -t 16
# Normal control sample
debreak --bam normal.bam -o output_normal -m min_supporting_reads --min_size 50 --ref reference_genome.fa -t 16
```

### 4\.4 SVIM

#### Germline SV Calling

```bash
svim alignment --minimum_depth min_supporting_reads --min_sv_size 30 output test.bam reference_genome.fa
bcftools view -i 'QUAL >= 5' variants.vcf > filtered_variants.vcf
```

#### Somatic SV Calling

```bash
# Tumor sample detection and filtering
svim alignment --minimum_depth min_supporting_reads --min_sv_size 50 output_tumor tumor.bam reference_genome.fa
bcftools view -i 'QUAL >= 5' tumor.vcf > filtered_tumor.vcf
# Normal control sample detection and filtering
svim alignment --minimum_depth min_supporting_reads --min_sv_size 50 output_normal normal.bam reference_genome.fa
bcftools view -i 'QUAL >= 5' normal.vcf > filtered_normal.vcf
```

### 4\.5 SVDF

#### Germline SV Calling

```bash
# CLR data
python svdf.py call test.bam --working_dir ./output -s min_supporting_reads --read_type CLR -t 16 --mode general --min_sv_size 30
# CCS data
python svdf.py call test.bam --working_dir ./output -s min_supporting_reads --read_type CCS -t 16 --mode general --min_sv_size 30
# ONT data
python svdf.py call test.bam --working_dir ./output -s min_supporting_reads --read_type ONT -t 16 --mode general --min_sv_size 30
```

#### Somatic SV Calling

```bash
# Tumor sample
python svdf.py call tumor.bam --working_dir ./output/tumor -s min_supporting_reads -t 16 --mode sensitive --min_sv_size 50
# Normal control sample
python svdf.py call normal.bam --working_dir ./output/normal -s min_supporting_reads -t 16 --mode sensitive --min_sv_size 50
```

### 4\.6 SVHunter

#### Germline SV Calling

```bash
python SVHunter.py generate ./test.bam ./datapath 16 includecontig
python SVHunter.py call ./predict_weight.h5 ./datapath ./test.bam ./predict_path ./outvcfpath 16 includecontig
```

#### Somatic SV Calling

```bash
# Tumor sample
python SVHunter.py generate ./tumor.bam ./datapath_tumor 16 includecontig
python SVHunter.py call ./predict_weight.h5 ./datapath_tumor ./tumor.bam ./predict_path ./outvcfpath 16 includecontig
# Normal control sample
python SVHunter.py generate ./normal.bam ./datapath_normal 16 includecontig
python SVHunter.py call ./predict_weight.h5 ./datapath_normal ./normal.bam ./predict_path ./outvcfpath 16 includecontig
```

### 4\.7 NanoVar

#### Germline SV Calling

```bash
# CLR data
nanovar -x pacbio-clr -c min_supporting_reads -t 16 -l 30 test.bam reference_genome.fa output
# CCS data
nanovar -x pacbio-ccs -c min_supporting_reads -t 16 -l 30 test.bam reference_genome.fa output
# ONT data
nanovar -x ont -c min_supporting_reads -t 16 -l 30 test.bam reference_genome.fa output
```

#### Somatic SV Calling

```bash
# CCS platform
nanovar -x pacbio-ccs -c min_supporting_reads -t 16 -l 50 tumor.bam reference_genome.fa tumor
nanovar -x pacbio-ccs -c min_supporting_reads -t 16 -l 50 normal.bam reference_genome.fa normal
# ONT platform
nanovar -x ont -c min_supporting_reads -t 16 -l 50 tumor.bam reference_genome.fa tumor
nanovar -x ont -c min_supporting_reads -t 16 -l 50 normal.bam reference_genome.fa normal
```

### 4\.8 SVision

#### Germline SV Calling

```bash
python SVision -o output -b test.bam -m ./svision_model/svision-cnn-model.ckpt -g reference_genome.fa -n test -s min_supporting_reads -t 16
```

#### Somatic SV Calling

```bash
# Tumor sample
python SVision -o tumor_output -b tumor.bam -m ./svision_model/svision-cnn-model.ckpt -g reference_genome.fa -n tumor -s min_supporting_reads -t 16
# Normal control sample
python SVision -o normal_output -b normal.bam -m ./svision_model/svision-cnn-model.ckpt -g reference_genome.fa -n normal -s min_supporting_reads -t 16
```

### 4\.9 SVision\-pro

#### Germline SV Calling

```bash
# CLR mode
python SVision-pro --target_path test.bam --genome_path reference_genome.fa --model_path ./src/pre_process/model_liteunet_256_8_16_32_32_32.pth --sample test --out_path output --detect_mode germline --process_num 16 --preset error-prone --min_supp min_supporting_reads
# CCS mode
python SVision-pro --target_path test.bam --genome_path reference_genome.fa --model_path ./src/pre_process/model_liteunet_256_8_16_32_32_32.pth --sample test --out_path output --detect_mode germline --process_num 16 --preset hifi --min_supp min_supporting_reads
# ONT mode
python SVision-pro --target_path test.bam --genome_path reference_genome.fa --model_path ./src/pre_process/model_liteunet_256_8_16_32_32_32.pth --sample test --out_path output --detect_mode germline --process_num 16 --preset error-prone --min_supp min_supporting_reads
```

#### Somatic SV Calling

```bash
# CCS mode
python SVision-pro --target_path tumor.bam --base_path normal.bam --model_path ./src/pre_process/model_liteunet_256_8_16_32_32_32.pth --out_path somatic_ouput --genome_path reference_genome.fa --sample_name somatic --preset hifi --process 16 --detect_mode somatic --min_supp min_supporting_reads
# ONT mode
python SVision-pro --target_path tumor.bam --base_path normal.bam --model_path ./src/pre_process/model_liteunet_256_8_16_32_32_32.pth --out_path somatic_ouput --genome_path reference_genome.fa --sample_name somatic --preset error-prone --process 16 --detect_mode somatic --min_supp min_supporting_reads
```

### 4\.10 Dysgu

#### Germline SV Calling

```bash
dysgu run -p 16 --min-size 30 --min-support min_supporting_reads reference_genome.fa output test.bam > test.vcf -x
```

#### Somatic SV Calling

```bash
# Tumor sample
dysgu run -p 16 --min-size 50 --min-support min_supporting_reads reference_genome.fa output tumor.bam > tumor.vcf -x
# Normal control sample
dysgu run -p 16 --min-size 50 --min-support min_supporting_reads reference_genome.fa output normal.bam > normal.vcf -x
```

### 4\.11 Delly2

#### Germline SV Calling

```bash
# CLR / CCS data
delly lr -t ALL -y pb -o test.bcf -g reference_genome.fa test.bam
# ONT data
delly lr -t ALL -y ont -o test.bcf -g reference_genome.fa test.bam
# Convert BCF to VCF
bcftools view test.bcf &gt; test.vcf
```

#### Somatic SV Calling

```bash
# CCS somatic calling
delly lr -t ALL -y pb -o somatic.bcf -g reference_genome.fa tumor.bam normal.bam
# ONT somatic calling
delly lr -t ALL -y ont -o somatic.bcf -g reference_genome.fa tumor.bam normal.bam
# Filter somatic variants
delly filter -f somatic -s sample.tsv -o somatic_filter.bcf somatic.bcf
bcftools view somatic_filter.bcf > somatic_filter.vcf
```

### 4\.12 SAVANA

#### Somatic SV Calling \(Tumor \+ Normal\)

```bash
savana --outdir output --tumour tumor.bam --normal normal.bam --ref reference_genome.fa --length 50 --threads 16
```

#### Tumor\-only SV Calling

```bash
savana to --tumour tumor.bam --outdir output --ref reference_genome.fa --g1000_vcf 1000g_hg38 --threads 16
```

### 4\.13 nanomonsv

#### Somatic SV Calling

```bash
# Parse tumor and normal bam files
nanomonsv parse tumor.bam output/test_tumor
nanomonsv parse normal.bam output/test_normal
# Detect somatic SVs with control background filtering
nanomonsv get output/test_tumor tumor.bam reference_genome.fa --control_prefix output/test_normal --control_bam normal.bam --use_racon --single_bnd --min_indel_size 50 --debug
```

### 4\.14 Severus

#### Germline SV Calling

```bash
python severus.py --target-bam test.bam --vntr-bed ./vntrs/human_hs37d5.bed --out-dir output -t 16 --min-support min_supporting_reads --min-sv-size 30
```

#### Somatic SV Calling \(Tumor \&amp; Normal Pair\)

```bash
python severus.py --target-bam tumor.bam --control-bam normal.bam --vntr-bed ./vntrs/human_GRCh38_no_alt_analysis_set.trf.bed.gz --out-dir output -t 16 --min-sv-size 50
```

#### Tumor\-only SV Calling

```bash
python severus.py --target-bam tumor.bam --out-dir output -t 16 --vntr-bed ./vntrs/human_GRCh38_no_alt_analysis_set.trf.bed.gz --PON ./pon/PoN_1000G_hg38.tsv.gz
```

---

## 5\. Somatic SV Filtering Strategies

### Tool Versions

jasmine \(v1\.1\.5\), bcftools

### 5\.1 Filter for Germline\-based SV Callers

```bash
# Merge multiple VCF outputs
ls *.vcf > filelist.txt
jasmine file_list=filelist.txt --dup_to_ins max_dist 500 genome_file=reference_genome.fa out_file=merged.vcf
# Extract pure somatic SVs (only present in tumor)
bcftools view -i "SUPP_VEC = '01'" merged.vcf > merge_somatic_only.vcf
```

### 5\.2 Filter for Sniffles2 Results

```bash
# Merge tumor and normal snf files
sniffles --input sniffles_normal.snf sniffles_tumor.snf --vcf merge_normal_tumor.vcf --allow-overwrite
# Extract somatic exclusive variations
bcftools view -i "SUPP_VEC = '01'" merge_normal_tumor.vcf > merge_somatic_only.vcf
```

---

## 6\. SV Evaluation

### Tool Versions

Truvari \(v5\.4\.0\), Minda \(v0\.0\.2\)

### 6\.1 Germline SV Evaluation \(Truvari\)

```bash
# Benchmark germline SV calls against gold standard dataset
truvari bench -f reference_genome.fa -b ./HG002_SVs_Tier1_v0.6.vcf.gz -o tool_eval --sizemin 50 --sizefilt 50 --passonly -p 0.00 -c variants.vcf.gz --includebed ./HG002_SVs_Tier1_v0.6.bed
```

### 6\.2 Somatic SV Evaluation \(Minda\)

```bash
# Evaluate somatic SV detection accuracy
python minda.py truthset --base truth.vcf --vcfs caller_vcf --min_size 50 --tolerance 500 --out_dir output
```

> （注：文档部分内容可能由 AI 生成）
