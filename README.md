# BrassiWood analysis pipeline

## Overview

This repository contains the scripts used to generate the analyses and figures presented in:

**Hendriks, K.P. *et al.* (2026)**
**Why woodiness repeatedly evolves—and disappears**
(submitted)

The workflow reproduces the principal analyses presented in the manuscript, including phylogenomic inference, ancestral-state reconstruction, environmental niche analyses, phylogenetic comparative analyses, and the generation of publication-quality figures.

---

## Repository structure

```text

BrassiWood-analysis/

├── WP1_BrassiToL/                                   Directory to create the Brassicaceae Tree of Life
    ├── WP1_BrassiToL/scripts/                       All scripts
    ├── WP1_BrassiToL/data/                          Required data
    ├── WP1_BrassiToL/results_intermediate/          Intermediate results
    ├── WP1_BrassiToL/results_final/                 Final results for publication
├── WP2_BrassiNiche/                                 Directory to analyse the growth form niches of woody and herbaceous Brassicaceae
    ├── WP1_BrassiToL/scripts/                       All scripts
    ├── WP1_BrassiToL/data/                          Required data
    ├── WP1_BrassiToL/results_intermediate/          Intermediate results
    ├── WP1_BrassiToL/results_final/                 Final results for publication
├── LICENSE
└── README.md

```



---

## Software requirements

The workflow relies on a combination of:

* Bash
* R
* Python
* HybPiper
* IQ-TREE
* ASTRAL
* BayesTraits
* BEAST
* TreeAnnotator

Additional package dependencies are documented within the individual scripts.

---

## Data availability

Large sequencing datasets are not included in this repository.

The datasets used in this study are available through the associated publication and the repositories cited therein.

---

## Citation

If you use this repository, please cite:

> Hendriks, K.P. *et al.* (2026).
> *Why woodiness repeatedly evolves—and disappears.*

Please also cite the archived Zenodo release associated with the version used.

---

## License

This repository is distributed under the MIT License. See the `LICENSE` file for details.

---

## Contact

**Kasper P. Hendriks**
Naturalis Biodiversity Center
Leiden, The Netherlands

Email: [k.p.hendriks@gmail.com](mailto:k.p.hendriks@gmail.com)
