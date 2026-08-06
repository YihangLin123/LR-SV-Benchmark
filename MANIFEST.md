# Formal Audit Manifest for LR-SV-Benchmark

| Benchmark Pipeline | Corresponding Download Links |
| :--- | :--- |
| **Dataset** | All datasets are available in Supplementary Table S2 |
| **Reference genome** | • GRCh37: http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz<br>• GRCh38: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa<br>• T2T-CHM13: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz |
| **Aligner** | • minimap2 (v2.30.0): https://github.com/lh3/minimap2<br>• pbmm2 (v26.1.99): https://github.com/PacificBiosciences/pbmm2<br>• ngmlr (v0.2.7): https://github.com/philres/ngmlr<br>*(Aligners for each dataset are available in Supplementary Table S4)* |
| **BAM file** | All BAM files are available in Supplementary Table S2 |
| **Caller version** | Detailed information for all callers, including versions and corresponding links, is available in Supplementary Table S1. |
| **Caller command** | Detailed running commands for all callers are available in Supplementary Note 6 and the GitHub repository at [https://github.com/model-lab/LR-SV-Benchmark](https://github.com/model-lab/LR-SV-Benchmark). |
| **Raw VCF** | The raw VCF files output by all callers are available in the GitHub repository at [https://github.com/model-lab/LR-SV-Benchmark](https://github.com/model-lab/LR-SV-Benchmark). |
| **Filtered VCF** | All filtered VCF files are available in the GitHub repository at [https://github.com/model-lab/LR-SV-Benchmark](https://github.com/model-lab/LR-SV-Benchmark). |
| **Ground truth** | All ground truth datasets are available in Supplementary Table S3. |
| **Benchmark tool and script** | • Truvari (v5.4.0): https://github.com/ACEnglish/truvari<br>• Minda (v0.0.2): https://github.com/KolmogorovLab/minda<br>*(All evaluation scripts are available in the GitHub repository)* |
| **Benchmark command** | Detailed running commands are available in Supplementary Note 6 and the GitHub repository. |
| **Final metric table** | All benchmarking results in this study are detailed in Supplementary Table S5-S69. |
