# Bird Microplastic Meta-Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## 📋 Overview

This repository contains data and code for a global meta-analysis investigating the ecological and evolutionary drivers of microplastic accumulation in birds. The study examines how morphological traits (body mass, beak length) and ecological factors (habitat, feeding guild, migration) predict interspecific variation in microplastic burden across 175 individuals representing 94 species.

**Key finding:** Smaller birds in terrestrial habitats carry the highest microplastic loads—a striking departure from the trophic transfer paradigm common in other taxa. A 164g body mass threshold and 31mm beak length cutoff consistently identify high-risk species.


## 🔬 Core Analyses

- **QA/QC Assessment**: Evaluation of methodological heterogeneity across extraction and identification methods
- **Phylogenetic Signal**: Pagel's λ for MP burden, body mass, and beak length
- **Decision Tree Analysis**: Identification of critical thresholds (164g, 31mm)
- **PGLS Models**: Trait-based predictors accounting for phylogeny
- **Functional Diversity**: Community-level vulnerability patterns
- **Tissue Accumulation**: Retention patterns across digestive organs

## 📊 Key Outputs

- Main text figures (1-5)
- Extended Data figures (ED1-ED6)
- Supplementary tables with full statistical outputs
- Phylogenetic tree and trait reconstructions

---
## 📦 Dependencies

Core R packages: tidyverse, ape, phytools, caper, ggtree, patchwork, viridis, rstatix
---