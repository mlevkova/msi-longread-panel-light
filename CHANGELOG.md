# Changelog

All notable changes to this project are documented in this file.

This project follows semantic versioning.

---

## [v1.0.0] - 2026-02-27

### Initial Stable Release

#### Added
- Targeted long-read transcript quantification pipeline
- Custom transcript panel alignment using minimap2
- RPM-based normalization workflow
- Multi-sample aggregation module
- Immune profiling module (T cell, cytotoxicity, checkpoint markers)
- Fusion-like multi-gene alignment detection
- Optional variant calling module (bcftools-based)
- Heatmap visualization example (Python)
- FFPE-compatible processing workflow
- Installation instructions with Conda environment
- Reproducibility documentation

#### Validated On
- FFPE colorectal carcinoma samples
- PCR-confirmed MSI-high tumors
- Oxford Nanopore MinION platform (R10.4.1 flow cell)
- cDNA library preparation

---

## Pre-release Development (2025–2026)

### Development Phase
- Initial MMR-only transcript panel
- Expansion to immune gene panel
- Optimization for low-resource workstations
- Debugging of barcode multiplex workflows
- Cross-run aggregation implementation
