# Bird Microplastic Meta-Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 📋 Overview

This repository contains data and code for a global meta-analysis investigating the ecological and evolutionary drivers of microplastic accumulation in birds. The study examines how morphological traits (body mass, wing length) and ecological factors (habitat, feeding guild, migration) predict interspecific variation in microplastic burden across 175 individuals representing 101 species.

## 🔬 Core Analyses

- **QA/QC Assessment**: Evaluation of methodological heterogeneity across extraction and identification methods
- **Phylogenetic Signal**: Pagel's λ for MP burden, body mass, and beak length
- **Decision Tree Analysis**: Identification of critical thresholds (164g, 31mm)
- **PGLS Models**: Trait-based predictors accounting for phylogeny
- **Functional Diversity**: Community-level vulnerability patterns
- **Tissue Accumulation**: Retention patterns across digestive organs

## 📊 Key Outputs

- Main text figures (1-3)
- Sensitivity analysis
- Supplementary tables with full statistical outputs
- Phylogenetic tree and trait reconstructions

---
## 📦 Dependencies

Core R packages: tidyverse, ape, phytools, caper, ggtree, patchwork, viridis, rstatix
---