# msi-longread-panel-light

Lightweight targeted long-read transcriptomic profiling pipeline for Oxford Nanopore cDNA sequencing data.

---

## Overview

This repository contains a reproducible bioinformatic workflow for targeted long-read transcriptomic analysis using Oxford Nanopore Technologies (ONT) sequencing data.

The pipeline is optimized for gene-panel–restricted expression profiling from cDNA libraries derived from formalin-fixed paraffin-embedded (FFPE) tissue, but can be applied to fresh-frozen samples as well.

The workflow enables:

- Targeted transcript-level quantification  
- Immune marker profiling  
- Mismatch repair (MMR) gene expression analysis  
- Lightweight normalization (RPM-based)  
- Multi-sample aggregation  
- Low computational footprint suitable for clinical environments  

The approach relies on targeted alignment to a custom transcript reference panel followed by direct count aggregation and normalization.

---

## Design Principles

This pipeline was developed with the following goals:

- Minimal computational requirements  
- Transparent and interpretable quantification  
- Compatibility with multiplexed barcode runs  
- Adaptability to custom transcript panels  
- Reproducibility in clinical research settings  

The workflow avoids heavy transcriptome-wide quantification frameworks and instead focuses on alignment-based targeted counting suitable for hypothesis-driven panels.

---

## Typical Use Case

This pipeline is suitable for:

- Targeted transcriptome analysis  
- Immune microenvironment profiling  
- Expression analysis of selected gene panels  
- FFPE-derived RNA sequencing  
- Small to medium-sized cohort studies  
- Clinical translational research projects  

Although originally designed for MSI-high colorectal cancer studies, it can be adapted to any transcript panel.

---

## Workflow Summary

1. Merge FASTQ files per barcode  
2. (Optional) Perform read orientation using PyChopper  
3. Align reads to custom transcript panel using `minimap2`  
4. Extract gene assignments from PAF output  
5. Generate raw gene counts  
6. Normalize counts to RPM (reads per million)  
7. Aggregate multi-sample results  
8. Perform immune-focused downstream analysis and visualization  

---

## Software Requirements

Required:

- Python ≥ 3.10  
- minimap2 ≥ 2.30  
- samtools  
- seqtk  
- awk  
- pandas  
- numpy  
- matplotlib  

Optional:

- bcftools (for exploratory variant analysis)
- pychopper (for ONT cDNA read orientation)

---

## Installation

We strongly recommend using Conda for environment management.

### 1. Create Environment

```bash
conda create -n ont-cdna python=3.12 -y
conda activate ont-cdna
```

### 2. Install Bioinformatics Tools

```bash
conda install -c bioconda minimap2=2.30 samtools seqtk bcftools pychopper -y
```

### 3. Install Python Packages

```bash
pip install pandas numpy matplotlib
```

### 4. Verify Installation

```bash
minimap2 --version
samtools --version
python --version
```

---

## Input Requirements

- Basecalled ONT FASTQ files (single or multiplexed)
- Custom transcript reference panel in FASTA format
- Barcoded directory structure (e.g., barcode01, barcode02, ...)

---

## Instructions (Usage)

Edit the two paths below before running:

```bash
set -euo pipefail
conda activate ont-cdna

BASE="/mnt/d/MSI final"
PANEL="/home/username/path/to/expanded_panel.labeled.fa"
```

### 1. Process All Barcode Folders

```bash
GENES='MLH1|MSH2|MSH6|PMS2|EPCAM|ACTB|GAPDH|RPLP0|RPS18|TP53|BRCA1|BRCA2|CD274|PDCD1|CTLA4'

OUT="$BASE/combined_counts_rpm.tsv"
echo -e "RunFolder\tSample\tGene\tCount\tRPM" > "$OUT"

for PASS in "$BASE"/fastq_pass*; do
  [ -d "$PASS" ] || continue
  runname=$(basename "$PASS")
  cd "$PASS"

  for B in barcode*; do
    [ -d "$B" ] || continue
    cd "$B"

    cat *.fastq *.fastq.gz 2>/dev/null > "${B}.merged.fastq" || true

    pychopper "${B}.merged.fastq" "${B}.pychop.fastq" || true

    minimap2 -K 64k -t 1 -x map-ont --secondary=no "$PANEL" "${B}.pychop.fastq" > "${B}.expanded.paf"

    awk -F'\t' '{split($6,a,"|"); g=a[1]; c[g]++}
    END{for(g in c) print g"\t"c[g]}' "${B}.expanded.paf" | sort > "${B}_counts.tsv"

    python3 - <<EOF
import pandas as pd
df=pd.read_csv("${B}_counts.tsv",sep="\t",header=None,names=["Gene","Count"])
tot=df["Count"].sum()
df["RPM"]=(df["Count"]*1e6/tot) if tot>0 else 0
df.to_csv("${B}_counts_rpm.tsv",sep="\t",index=False)
EOF

    awk -v run="$runname" -v samp="$B" -v pat="$GENES" 'BEGIN{FS=OFS="\t"}
    NR>1 && $1 ~ "^(" pat ")$" {print run,samp,$1,$2,$3}' "${B}_counts_rpm.tsv" >> "$OUT"

    cd ..
  done
done

echo "Pipeline completed successfully."
```

---

## Output Files

For each barcode:

- `*_expanded.paf` – alignment file  
- `*_counts.tsv` – raw gene counts  
- `*_counts_rpm.tsv` – normalized expression table  

Global output:

- `combined_counts_rpm.tsv` – aggregated multi-sample expression matrix  

---

## Reproducibility

All analyses were performed using:

- Python 3.12  
- minimap2 v2.30  
- Single-threaded alignment mode  
- RPM-based normalization  

---

## Version

Current stable release: v1.0.0

---

## Citation

If you use this pipeline, please cite the associated publication (details to be updated upon acceptance). Please contact Mariya Levkova.

---

## License

MIT License
