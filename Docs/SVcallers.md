# Structural Variant callers

| SV callers | Version | SV types | Target | Release year |
| :--- | :---: | :--- | :---: | :---: |
| **[Delly2](https://github.com/dellytools/delly)** | 1.5.0 | DEL, INS, DUP, INV, TRA | Both | 2012 |
| **[SVIM](https://github.com/eldariont/svim)** | 2.0.0 | DEL, INS, DUP, INV, BND, DUP:TANDEM, DUP:INT | Germline | 2019 |
| **[NanoVar](https://github.com/cytham/nanovar)** `[DL]` | 1.8.3 | DEL, INS, DUP, INV, BND | Germline | 2020 |
| **[cuteSV2](https://github.com/tjiangHIT/cuteSV)** | 2.1.3 | DEL, INS, DUP, INV, BND | Germline | 2020 |
| **[Dysgu](https://github.com/kcleal/dysgu)** `[ML]` | 1.7.0 | DEL, INS, DUP, INV, TRA | Germline | 2022 |
| **[SVision](https://github.com/xjtu-omics/SVision)** `[DL]` | 1.4.0 | SSVs (DEL, INS, DUP, INV, DUP:TANDEM)<br>CSVs (25 types, e.g., DEL+INV, INV+DUP) | Germline | 2022 |
| **[DeBreak](https://github.com/Maggi-Chen/DeBreak)** | 1.0.2 | DEL, INS, DUP, INV, TRA | Germline | 2023 |
| **[nanomonsv](https://github.com/friend1ws/nanomonsv)** | 0.8.0 | DEL, INS, DUP, INV, BND | Somatic | 2023 |
| **[Sniffles2](https://github.com/fritzsedlazeck/Sniffles)** | 2.6.3 | DEL, INS, DUP, INV, BND | Germline | 2024 |
| **[SVision-pro](https://github.com/songbowang125/SVision-pro)** `[DL]` | 2.5.0 | SSVs (DEL, INS, DUP, INV, BND, DUP:TANDEM)<br>CSVs (Multi-component combinations) | Both | 2024 |
| **[SVDF](https://github.com/coopsor/SVDF)** `[DL]` | default | DEL, INS, DUP, INV, TRA | Germline | 2024 |
| **[SVHunter](https://github.com/eioyuou/SVHunter)** `[DL]` | 1.1.0 | DEL, INS, DUP, INV, TRA | Germline | 2025 |
| **[SAVANA](https://github.com/cortes-ciriano-lab/savana)** | 1.6.0 | DEL, INS, DUP, INV, BND, SCNA | Somatic | 2025 |
| **[Severus](https://github.com/KolmogorovLab/Severus)** | 1.6.0 | DEL, INS, DUP, INV, BND | Both | 2026 |

---

## Acronyms & Abbreviations

* **SV Types:**
  * `DEL`: Deletion
  * `INS`: Insertion
  * `DUP`: Duplication (`DUP:TANDEM`: Tandem Duplication, `DUP:INT`: Interspersed Duplication)
  * `INV`: Inversion
  * `TRA`: Translocation
  * `BND`: Breakend / Complex Rearrangement
  * `SCNA`: Somatic Copy Number Alteration
  * `SSVs`: Simple Structural Variants
  * `CSVs`: Complex Structural Variants
* **Methodology Badges:**
  * `[ML]`: Machine Learning based algorithm
  * `[DL]`: Deep Learning based model
