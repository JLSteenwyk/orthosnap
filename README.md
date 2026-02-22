<p align="center">
  <a href="https://github.com/jlsteenwyk/orthosnap">
    <img src="https://raw.githubusercontent.com/JLSteenwyk/orthosnap/master/docs/_static/img/logo.jpg" alt="Logo" width="400">
  </a>
  <p align="center">
    <a href="https://jlsteenwyk.com/orthosnap/">Docs</a>
    ·
    <a href="https://github.com/jlsteenwyk/orthosnap/issues">Report Bug</a>
    ·
    <a href="https://github.com/jlsteenwyk/orthosnap/issues">Request Feature</a>
  </p>
    <p align="center">
        <a href="https://github.com/JLSteenwyk/orthosnap/actions" alt="Build">
            <img src="https://img.shields.io/github/actions/workflow/status/JLSteenwyk/orthosnap/ci.yml?branch=master">
        </a>
        <a href="https://codecov.io/gh/JLSteenwyk/orthosnap">
          <img src="https://codecov.io/gh/JLSteenwyk/orthosnap/branch/master/graph/badge.svg?token=FX66FUET0L"/>
        </a>
        <a href="https://github.com/JLSteenwyk/orthosnap/graphs/contributors" alt="Contributors">
            <img src="https://img.shields.io/github/contributors/JLSteenwyk/orthosnap">
        </a>
        <a href="https://bsky.app/profile/jlsteenwyk.bsky.social" target="_blank" rel="noopener noreferrer">
          <img src="https://img.shields.io/badge/Bluesky-0285FF?logo=bluesky&logoColor=fff">
        </a>
        <br />
        <a href="https://pepy.tech/badge/orthosnap">
          <img src="https://static.pepy.tech/personalized-badge/orthosnap?period=total&units=international_system&left_color=grey&right_color=blue&left_text=PyPi%20Downloads">
        </a>
        <a href="https://lbesson.mit-license.org/" alt="License">
            <img src="https://img.shields.io/badge/License-MIT-blue.svg">
        </a>
        <a href="https://pypi.org/project/orthosnap/" alt="PyPI - Python Version">
            <img src="https://img.shields.io/pypi/pyversions/orthosnap">
        </a>
        <a href="https://jlsteenwyk.com/publication_pdfs/2022_Steenwyk_etal_PLoS_Biology.pdf">
          <img src="https://zenodo.org/badge/DOI/10.1371/journal.pbio.3001827.svg">
        </a>
    </p>
</p>

OrthoSNAP is a tree splitting and pruning tool for retrieving single-copy orthologous subgroups (SNAP-OGs) from larger gene families.

If you found OrthoSNAP useful, please cite:
*OrthoSNAP: a tree splitting and pruning algorithm for retrieving single-copy orthologs from gene family trees*. Steenwyk et al. 2022, PLOS Biology. DOI: [10.1371/journal.pbio.3001827](https://jlsteenwyk.com/publication_pdfs/2022_Steenwyk_etal_PLoS_Biology.pdf).

---

Full usage documentation and tutorial:
[https://jlsteenwyk.com/orthosnap/](https://jlsteenwyk.com/orthosnap/)

## What's new in v1.6.0

Compared to v1.5.0 (plotting + performance improvements), v1.6.0 adds workflow-scale and reproducibility features:

- `--manifest`: batch execution from TSV/CSV manifests.
- `--validate-only`: preflight input concordance checks without extraction.
- `--structured-output`: machine-readable run metadata (`.run.json`) and subgroup summaries (`.subgroups.tsv`).
- `--occupancy-count` / `--occupancy-fraction`: explicit occupancy semantics.
- `--resume`: skip rerunning completed analyses.
- `--bootstrap-trees` + `--consensus-min-frequency` + `--consensus-trees`: consensus subgrouping across bootstrap tree uncertainty.

Compared to older releases:

- v1.5.0 focused on plotting and runtime optimization.
- v1.3.2 introduced configurable delimiters.
- v1.2.0 added inparalog handling reports.
- v1.0.0 and earlier focused on core pruning behavior.

## Installation

### Install with pip (recommended)

```shell
python -m venv .venv
source .venv/bin/activate
pip install orthosnap
```

### Install from source

```shell
git clone https://github.com/JLSteenwyk/orthosnap.git
cd orthosnap
python -m venv .venv
source .venv/bin/activate
make install
```

### Install with conda

```shell
conda install -c jlsteenwyk orthosnap
```

Conda package details:
https://anaconda.org/jlsteenwyk/orthosnap

## Quick start

```shell
orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre
```

Generate a color-coded SNAP-OG assignment plot for the full tree:

```shell
orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre -ps
```

Choose plot format (`png` default, `pdf` or `svg`):

```shell
orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre -ps -pf svg
```

Show all CLI options:

```shell
orthosnap -h
```

Run validation checks only (no subgroup extraction):

```shell
orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre --validate-only
```

Write structured provenance outputs (`.run.json` and `.subgroups.tsv`):

```shell
orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre --structured-output
```

Resume an interrupted or previously completed run:

```shell
orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre --resume
```

Use explicit occupancy semantics:

```shell
orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre --occupancy-count 5
orthosnap -f orthogroup_of_genes.faa -t phylogeny_of_orthogroup_of_genes.tre --occupancy-fraction 0.5
```

Run many orthogroups from a manifest (TSV/CSV with `tree` and `fasta` columns):

```shell
orthosnap --manifest runs.tsv --structured-output -op results/
```

Run bootstrap consensus mode using a file of tree paths (one per line):

```shell
orthosnap -f orthogroup_of_genes.faa -t reference.treefile --bootstrap-trees bootstrap_paths.txt --consensus-min-frequency 0.5
```

Also write consensus Newick trees:

```shell
orthosnap -f orthogroup_of_genes.faa -t reference.treefile --bootstrap-trees bootstrap_paths.txt --consensus-trees
```

## Support

If installation fails in a clean virtual environment, contact Jacob L. Steenwyk via:
- Email: https://jlsteenwyk.com/contact.html
- Twitter/X: https://twitter.com/jlsteenwyk
