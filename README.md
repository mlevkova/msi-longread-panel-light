# msi-longread-panel-light
Lightweight targeted long-read transcriptomic profiling pipeline for Oxford Nanopore cDNA sequencing data (FFPE-optimized).

---

## Overview

This repository contains a reproducible bioinformatic workflow for targeted long-read transcriptomic analysis using Oxford Nanopore Technologies (ONT) sequencing data.  
It is designed for gene-panel–restricted expression profiling from cDNA libraries obtained from **formalin-fixed paraffin-embedded (FFPE)** tumor tissue.

The workflow supports:

- Targeted transcript quantification  
- Mismatch repair (MMR) gene profiling  
- Immune marker profiling  
- Detection of fusion-like multi-gene alignments  
- Lightweight normalization (RPM)  
- Cohort-wide aggregation  
- Optional variant calling for exploratory analysis  

The pipeline is computationally lightweight and runs on standard laboratory workstations without HPC resources.

---

## Features

### Targeted transcript quantification  
Fast alignment against a curated transcript panel using minimap2.

### Immune microenvironment scoring  
Detection of:
- T cell activation genes  
- Cytotoxicity genes (NKG7, PRF1, GZMB, IFNG)  
- Immune checkpoints (PDCD1, CD274, CTLA4, HAVCR2, LAG3, TIGIT)  

### MMR profiling  
Quantification of:
- MLH1  
- MSH2  
- MSH6  
- PMS2  
- EPCAM  

### Fusion-like read detection  
Identification of reads aligning to *multiple* target genes, supporting:
- Complex rearrangements  
- Potential fusion transcripts  
- Artifacts in highly rearranged MSI tumors  

### Visualization module  
Optional Python scripts generate:
- Heatmaps  
- Barplots  
- Immune activation profiles  
- Quality control metrics  

### FFPE compatibility  
Validated on multiple runs using FFPE colorectal carcinoma tissue.

---

## Typical Use Cases

This pipeline is ideal for:

- Targeted long-read transcriptomics  
- Clinical translational research  
- Immune microenvironment characterization  
- FFPE sequencing projects  
- MSI-high colorectal cancer studies  
- Low-input or degraded RNA workflows  

---

## Workflow Summary

1. Merge FASTQ files for each barcode  
2. (Optional) Process with PyChopper to orient cDNA reads  
3. Align reads to custom transcript panel  
4. Extract gene names from minimap2 PAF  
5. Generate raw counts per gene  
6. Normalize to RPM  
7. Aggregate across multiple sequencing runs  
8. Perform one or more optional analyses:
   - Immune profiling  
   - Fusion read detection  
   - Variant calling  
   - Heatmap visualization  

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

- **pychopper** (for read orientation)
- **bcftools** (for variant calling)

---

## Installation

### 1. Create a conda environment
```bash
conda create -n ont-cdna python=3.12 -y
conda activate ont-cdna
```

### 2. Install core tools
```bash
conda install -c bioconda minimap2=2.30 samtools seqtk bcftools pychopper -y
```

### 3. Install Python dependencies
```bash
pip install pandas numpy matplotlib seaborn
```

---

## Input Requirements

- Basecalled ONT FASTQ files  
- Barcoded directory structure:  
  `fastq_pass/run/barcode01/`  
- Custom transcript reference (FASTA)  
  - Example included:  
    `expanded_panel.labeled.fa`

---

# Main Pipeline (Expression Quantification)

Edit these two paths first:

```bash
BASE="/mnt/d/MSI final"
PANEL="/home/username/path/to/expanded_panel.labeled.fa"
conda activate ont-cdna
set -euo pipefail
```

Genes collected:

```bash
GENES='MLH1|MSH2|MSH6|PMS2|EPCAM|ACTB|GAPDH|RPLP0|RPS18|TP53|BRCA1|BRCA2|CD274|PDCD1|CTLA4'
```

### Run the pipeline
```bash
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

    minimap2 -K 64k -t 1 -x map-ont --secondary=no "$PANEL" \
      "${B}.pychop.fastq" > "${B}.expanded.paf"

    awk -F'\t' '{split($6,a,"|"); g=a[1]; c[g]++}
    END{for(g in c) print g"\t"c[g]}' "${B}.expanded.paf" \
      | sort > "${B}_counts.tsv"

    python3 - <<EOF
import pandas as pd
df=pd.read_csv("${B}_counts.tsv",sep="\t",header=None,names=["Gene","Count"])
tot=df["Count"].sum()
df["RPM"]=(df["Count"]*1e6/tot) if tot>0 else 0
df.to_csv("${B}_counts_rpm.tsv",sep="\t",index=False)
EOF

    awk -v run="$runname" -v samp="$B" -v pat="$GENES" \
      'BEGIN{FS=OFS="\t"} NR>1 && $1 ~ "^(" pat ")$" {print run,samp,$1,$2,$3}' \
      "${B}_counts_rpm.tsv" >> "$OUT"

    cd ..
  done
done
```

---

# Optional Analysis Modules

## 1. Immune Profiling Module
Extract:

- T cell markers (CD3D, CD3E, CD8A, TRAC)  
- Cytotoxic markers (NKG7, PRF1, GZMB, IFNG)  
- Immune checkpoints (PDCD1, CD274, CTLA4, HAVCR2)

```bash
awk -F'\t' '$3~/^(CD3D|CD3E|CD8A|TRAC|NKG7|PRF1|GZMB|IFNG|PDCD1|CD274|CTLA4|HAVCR2)$/' \
combined_counts_rpm.tsv | column -t
```

---

## 2. Fusion-like Read Detection
```bash
awk -F'\t' '{rid=$1; split($6,a,"|"); gene=a[1]; print rid"\t"gene}' \
  *.expanded.paf | sort -u > read_to_gene.tsv

cut -f1,2 read_to_gene.tsv | sort -u \
 | awk '{k[$1][$2]=1} END{for(r in k){c=0; g=""; for (x in k[r]){c++; g=g","x} if (c>1) print r"\t"c"\t"substr(g,2)}}' \
 > fusion_candidates.tsv
```

---

## 3. Variant Calling (Exploratory)
```bash
minimap2 -ax splice -t 2 GRCh38.fa sample.fastq | samtools sort -o aln.bam
samtools index aln.bam

bcftools mpileup -f GRCh38.fa aln.bam | bcftools call -mv -Oz -o variants.vcf.gz
```

---

## 4. Heatmap Visualization (Python)

Example:

```python
import pandas as pd, seaborn as sns, matplotlib.pyplot as plt

df = pd.read_csv("combined_counts_rpm.tsv", sep="\t")
pivot = df.pivot_table(index="Gene", columns="Sample", values="RPM")

sns.clustermap(pivot, cmap="viridis", figsize=(8,12))
plt.savefig("heatmap.png", dpi=300)
```

---

## Output Files

Per barcode:
- `*.expanded.paf`
- `*_counts.tsv`
- `*_counts_rpm.tsv`

Global:
- `combined_counts_rpm.tsv`
- Optional: heatmaps, fusion tables, variant VCF

---

## Reproducibility

- Python 3.12  
- minimap2 2.30  
- Single-threaded execution  
- RPM normalization  

---

## Version

Stable release: **v1.0.0**

---

## Citation

Please cite this repository when publishing results.  
Manuscript details will be added after acceptance.

---

## DOI

Archived at Zenodo:
https://doi.org/10.5281/zenodo.18824372

---

## License

MIT License
