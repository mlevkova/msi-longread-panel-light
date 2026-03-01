# msi-longread-panel-light

Lightweight targeted long-read transcriptomic profiling pipeline for Oxford Nanopore cDNA data.

---

## Overview

This repository contains a reproducible bioinformatic workflow for targeted long-read transcriptomic analysis using Oxford Nanopore Technologies (ONT) sequencing data.

The pipeline is optimized for gene-panel–restricted expression profiling from cDNA libraries, including challenging clinical material such as formalin-fixed paraffin-embedded (FFPE) tissue.

It provides:

- Targeted transcript-level quantification  
- Immune marker profiling  
- Mismatch repair (MMR) gene expression analysis  
- Lightweight normalization (RPM-based)  
- Low computational footprint suitable for clinical environments  

The workflow is alignment-based and designed to run on standard laboratory workstations without high-performance computing infrastructure.

---

## Design Principles

This pipeline was developed with the following goals:

- Minimal computational requirements  
- Transparent and interpretable quantification  
- Compatibility with multiplexed barcode runs  
- Adaptability to custom transcript panels  
- Reproducibility in clinical research settings  

The approach relies on targeted alignment to a custom transcript reference panel followed by direct count aggregation and normalization.

---

## Typical Use Case

The workflow is particularly suited for:

- Targeted transcriptome analysis  
- Immune microenvironment profiling  
- Expression analysis of selected gene panels  
- FFPE-derived RNA sequencing  
- Small to medium-sized cohort studies  

Although originally designed for mismatch repair–deficient colorectal cancer research, the pipeline can be adapted to any custom transcript panel.

---

## Workflow Summary

1. Merge FASTQ files per barcode  
2. Align reads to custom transcript panel using `minimap2`  
3. Extract gene assignments from PAF output  
4. Generate raw gene counts  
5. Normalize counts to RPM (reads per million)  
6. Aggregate multi-sample results  
7. Optional immune-focused analysis and visualization  

---

## Software Requirements

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

---

## Installation

Create a conda environment:

```bash
conda create -n ont-cdna python=3.12 -y
conda activate ont-cdna
conda install -c bioconda minimap2=2.30 samtools seqtk bcftools -y
pip install pandas numpy matplotlib
```

Verify installation:

```bash
minimap2 --version
samtools --version
python --version
```

---

## Instructions (Usage)

Edit the two paths below before running.

```bash
set -euo pipefail
conda activate ont-cdna

# =========================
BASE="/mnt/d/MSI final"
PANEL="/home/maria_levkova/barcode01/mmr_paf/expanded_panel.labeled.fa"
# =========================

GENES='MLH1|MSH2|MSH6|PMS2|EPCAM|ACTB|GAPDH|RPLP0|RPS18|TP53|BRCA1|BRCA2|CD274|PDCD1|CTLA4'

OUT="$BASE/combined_counts_rpm.tsv"
echo -e "RunFolder\tSample\tGene\tCount\tRPM" > "$OUT"

[ -s "$PANEL" ] || { echo "Panel FASTA not found"; exit 1; }

for PASS in "$BASE"/fastq_pass*; do
  [ -d "$PASS" ] || continue
  runname=$(basename "$PASS")
  cd "$PASS"

  for B in barcode*; do
    [ -d "$B" ] || continue

    sz_mb=$(du -sm "$B" | awk '{print $1}')
    if [ "$sz_mb" -lt 5 ]; then
      continue
    fi

    cd "$B"

    MERGED="${B}.all.fastq"
    : > "$MERGED"
    for f in *.fastq *.fq *.fastq.gz *.fq.gz 2>/dev/null; do
      [ -e "$f" ] || continue
      if [[ "$f" == *.gz ]]; then
        zcat "$f" >> "$MERGED"
      else
        cat "$f" >> "$MERGED"
      fi
    done

    PCH="${B}.pychop.full.fastq"
    pychopper "$MERGED" "$PCH"

    PAF="${B}.expanded.paf"
    minimap2 -K 64k -t 1 -x map-ont --secondary=no "$PANEL" "$PCH" > "$PAF"

    COUNTS="${B}_counts.tsv"
    awk -F'\t' '{split($6,a,"|"); g=a[1]; c[g]++}
    END{for(g in c) print g"\t"c[g]}' "$PAF" | sort > "$COUNTS"

    RPM="${B}_counts_rpm.tsv"
    python3 - <<EOF
import pandas as pd
df=pd.read_csv("$COUNTS",sep="\t",header=None,names=["Gene","Count"])
tot=df["Count"].sum()
df["RPM"]=(df["Count"]*1e6/tot) if tot>0 else 0
df.sort_values("Gene").to_csv("$RPM",sep="\t",index=False)
EOF

    awk -v run="$runname" -v samp="$B" -v pat="$GENES" 'BEGIN{FS=OFS="\t"}
      NR>1 && $1 ~ "^(" pat ")$" {print run,samp,$1,$2,$3}' "$RPM" >> "$OUT"

    cd ..
  done
done

echo "DONE"
column -t "$OUT" | head
```

---

## License

MIT License
