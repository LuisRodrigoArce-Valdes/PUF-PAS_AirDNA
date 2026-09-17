## Welcome to the repository for the bioinformatic analyses of the *Toward Repurposing Global Passive Air SamplingNetworks for Insect Monitoring: Promises and Pitfalls of Airborne eDNA* manuscript by Vilanova et al. 2026.

### Repository contents
This project analysed 2 libraries totaling 59 samples of airborne eDNA sampled with polyurethane foam passive air samplers. A COI fragment was amplified using the Leray set of primers (Morrill et al. 2021) and sequenced on a Miseq using a MiSeq reagent kit v3 (600-cycles; Illumina). This repository includes the scripts to perform all the bioinformatic processing of these samples.

### Repository structure:

This repository splits the main steps of the bioinformatic processing into multiple directories. Within each directory, subdirectories may be found, which include main steps across the QIIME2 pipeline. All directories and subdirectories are ordered in the way they should be run. Finally, each subdirectory branches into `bin` directories that contain the scripts associated to each step.

The `00_Databases` and `00_RawData` directories are not tracked by GitHub because of their size. They contain the raw data used for the analyses. Data should be downloaded and included in a structure as is shown in the repository tree to follow the pipeline.

### Repository tree:

```
.
├── 00_COI_Database                                 *Scrips to process the COI reference database.
│   └── bin
├── 00_Databases                                    *This directory contains and procceses the [COInr](https://zenodo.org/records/15515860) database that we used for taxonomic asignation*
│   ├── COInr_2025_05_23
│   ├── COInrLeray
│   └── COInrLerayQIIME
├── 00_RawData                                      *This directory contains the raw sequences split between the two sequencing libraries.
│   ├── PUFs_Pt1
│   └── PUFs_Pt2
├── 01_QIIME_pt1                                    *Scripts and intermediary steps for the processing of the first library.
│   ├── 01_Deblur
│   │   ├── bin
│   │   │   └── logs
│   │   ├── data
│   │   └── results
│   │       └── 00_FastQC
│   ├── 02_DADA2
│   │   ├── bin
│   │   │   └── logs
│   │   ├── data
│   │   └── results
│   └── 03_NaiveBayes
│       ├── bin
│       ├── data
│       └── results
│           └── 03_Matrices
│               ├── DADA2_t0_l10
│               │   ├── css
│               │   ├── js
│               │   └── q2templateassets
│               │       ├── css
│               │       ├── fonts
│               │       ├── img
│               │       └── js
│               ├── DADA2_t20_l10
│               │   ├── css
│               │   ├── js
│               │   └── q2templateassets
│               │       ├── css
│               │       ├── fonts
│               │       ├── img
│               │       └── js
│               ├── DADA2_t20_l50
│               │   ├── css
│               │   ├── js
│               │   └── q2templateassets
│               │       ├── css
│               │       ├── fonts
│               │       ├── img
│               │       └── js
│               ├── DADA2_t26_l10
│               │   ├── css
│               │   ├── js
│               │   └── q2templateassets
│               │       ├── css
│               │       ├── fonts
│               │       ├── img
│               │       └── js
│               ├── Deblur_t0_l10
│               │   ├── css
│               │   ├── js
│               │   └── q2templateassets
│               │       ├── css
│               │       ├── fonts
│               │       ├── img
│               │       └── js
│               ├── Deblur_t20_l10
│               │   ├── css
│               │   ├── js
│               │   └── q2templateassets
│               │       ├── css
│               │       ├── fonts
│               │       ├── img
│               │       └── js
│               ├── Deblur_t20_l50
│               │   ├── css
│               │   ├── js
│               │   └── q2templateassets
│               │       ├── css
│               │       ├── fonts
│               │       ├── img
│               │       └── js
│               └── Deblur_t26_l10
│                   ├── css
│                   ├── js
│                   └── q2templateassets
│                       ├── css
│                       ├── fonts
│                       ├── img
│                       └── js
├── 01_QIIME_pt2                                     *Scripts and intermediary steps for the processing of the second library.
│   ├── 01_Deblur
│   │   ├── bin
│   │   │   └── logs
│   │   ├── data
│   │   └── results
│   │       └── 00_FastQC
│   ├── 02_DADA2
│   │   ├── bin
│   │   │   └── logs
│   │   ├── data
│   │   └── results
│   └── 03_NaiveBayes
│       ├── bin
│       ├── data
│       └── results
│           └── 03_Matrices
│               ├── DADA2_t0_l10
│               │   ├── css
│               │   ├── js
│               │   └── q2templateassets
│               │       ├── css
│               │       ├── fonts
│               │       ├── img
│               │       └── js
│               ├── DADA2_t20_l10
│               │   ├── css
│               │   ├── js
│               │   └── q2templateassets
│               │       ├── css
│               │       ├── fonts
│               │       ├── img
│               │       └── js
│               ├── DADA2_t20_l50
│               │   ├── css
│               │   ├── js
│               │   └── q2templateassets
│               │       ├── css
│               │       ├── fonts
│               │       ├── img
│               │       └── js
│               ├── DADA2_t26_l10
│               │   ├── css
│               │   ├── js
│               │   └── q2templateassets
│               │       ├── css
│               │       ├── fonts
│               │       ├── img
│               │       └── js
│               ├── Deblur_t0_l10
│               │   ├── css
│               │   ├── js
│               │   └── q2templateassets
│               │       ├── css
│               │       ├── fonts
│               │       ├── img
│               │       └── js
│               ├── Deblur_t20_l10
│               │   ├── css
│               │   ├── js
│               │   └── q2templateassets
│               │       ├── css
│               │       ├── fonts
│               │       ├── img
│               │       └── js
│               ├── Deblur_t20_l50
│               │   ├── css
│               │   ├── js
│               │   └── q2templateassets
│               │       ├── css
│               │       ├── fonts
│               │       ├── img
│               │       └── js
│               └── Deblur_t26_l10
│                   ├── css
│                   ├── js
│                   └── q2templateassets
│                       ├── css
│                       ├── fonts
│                       ├── img
│                       └── js
├── 02_Filtering                                        *ASVs filtering schemes directory.
│   ├── bin
│   ├── data
│   └── results
└── meta

```

### Additional analyses:
The raw data can be downloaded from the NCBI accession number: [PRJNA1436621](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1436621).
This pipeline as well as the scripts of the implemented statistical analyses can be found on [figshare](https://figshare.com/s/01a97cef04a571adfeea).

