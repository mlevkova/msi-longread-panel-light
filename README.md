# msi-longread-panel-light
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
conda create -n ont-cdna python=3.12
conda activate ont-cdna
conda install -c bioconda minimap2 samtools seqtk
pip install pandas numpy matplotlib
