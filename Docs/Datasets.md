# Datasets

## Contents
1. Germline
2. Somatic
3. Benchmark truthsets

---

## 1. Germline

| Sample | Technology | Coverage | Note / Source |
| :--- | :--- | :--- | :--- |
| **HG002** | [CLR](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/) | 69× | Ashkenazim Trio (Son) |
| | [CCS](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb_20kb_chemistry2/) | 28× | |
| | [ONT](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/) | 48× | |
| **HG003** | [CLR](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/PacBio_MtSinai_NIST/) | 30× | Ashkenazim Trio (Father) |
| | [CCS](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/PacBio_CCS_15kb_20kb_chemistry2/) | 30× | |
| | [ONT](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/UCSC_Ultralong_OxfordNanopore_Promethion/) | 30× | |
| **HG004** | [CLR](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/PacBio_MtSinai_NIST/) | 30× | Ashkenazim Trio (Mother) |
| | [CCS](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/PacBio_CCS_15kb_20kb_chemistry2/) | 30× | |
| | [ONT](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/UCSC_Ultralong_OxfordNanopore_Promethion/) | 30× | |
| **HG005** | [CLR](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG005_NA24631_son/MtSinai_PacBio/) | 30× | Chinese Trio (Son) |
| | [CCS](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG005_NA24631_son/PacBio_CCS_15kb_20kb_chemistry2/) | 30× | |
| | [ONT](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG005_NA24631_son/UCSC_Ultralong_OxfordNanopore_Promethion/) | 30× | |
| **HG006** | [CLR](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG006_NA24694-huCA017E_father/PacBio_MtSinai/) | 30× | Chinese Trio (Father) |
| | [CCS](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG006_NA24694-huCA017E_father/PacBio_CCS_15kb_20kb_chemistry2/) | 30× | |
| | [ONT](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG006_NA24694-huCA017E_father/UCSC_Ultralong_OxfordNanopore_Promethion/) | 30× | |
| **HG007** | [CLR](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG007_NA24695-hu38168_mother/PacBio_MtSinai/) | 30× | Chinese Trio (Mother) |
| | [CCS](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG007_NA24695-hu38168_mother/PacBio_CCS_15kb_20kb_chemistry2/) | 30× | |
| | [ONT](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/ChineseTrio/HG007_NA24695-hu38168_mother/UCSC_Ultralong_OxfordNanopore_Promethion/) | 30× | |
| **CHM13** | [CLR](https://github.com/marbl/CHM13) | 54× | T2T Consortium |
| | [CCS](https://github.com/marbl/CHM13) | 33× | |
| | [ONT](https://github.com/marbl/CHM13) | 123× | |
| **SIM** | CLR | 30× | Simulated from CHM1 (nstd137) & KWS1 (nstd106) via [VISOR](https://github.com/davidebolo1993/VISOR) |
| | CCS | 30× | |
| | ONT | 30× | |

## 2. Somatic

| Sample | Type | Technology | Coverage | Source |
| :--- | :--- | :--- | :--- | :--- |
| **HCC1395** | Tumor | [CCS](https://downloads.pacbcloud.com/public/revio/2023Q2/HCC1395/HCC1395/) | 62× | PacBio Revio |
| | | [ONT](http://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG) | 19× | SEQC2 |
| | Normal | [CCS](https://downloads.pacbcloud.com/public/revio/2023Q2/HCC1395/HCC1395-BL/) | 43× | PacBio Revio |
| | | [ONT](http://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG) | 12× | SEQC2 |
| **COLO829** | Tumor | [CCS](https://downloads.pacbcloud.com/public/revio/2024Q4/WGS/COLO829/COLO829/) | 66× | PacBio Revio |
| | | [ONT](https://epi2me.nanoporetech.com/colo-2024.03) | 49× | ONT EPI2ME |
| | Normal | [CCS](https://downloads.pacbcloud.com/public/revio/2024Q4/WGS/COLO829/COLO829BL/) | 67× | PacBio Revio |
| | | [ONT](https://epi2me.nanoporetech.com/colo-2024.03) | 46× | ONT EPI2ME |
| **HG008** | Tumor | [CCS](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/PacBio_Revio_20240125/) | 84× | GIAB Somatic |
| | | [ONT](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/UCSC_ONT_20231003/) | 47× | GIAB Somatic |
| | Normal | [CCS](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/PacBio_Revio_20240125/) | 36× | GIAB Somatic |
| | | [ONT](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/Northeastern_ONT-std_20240422/) | 41× | GIAB Somatic |
| **Synthetic** | Tumor | [CCS](https://www.ncbi.nlm.nih.gov/sra/SRX10759866) | 63× | CHM1 / [CHM13](https://github.com/marbl/CHM13) Synthetic |
| | Normal | CCS | 30× | CHM1 |

## 3. Benchmark truthsets
