# Computergestuetzte Strukturelle Biologie

This repository contains a molecular dynamics study of the human orexin-1 receptor (OX1R) in complex with the selective antagonist GSK1059865. It includes the prepared simulation system, GROMACS input files, analysis scripts, selected derived results, and the Quarto sources for the final written report and presentation material.

## Project Focus

The main scientific goal of the project is to characterize the structural stability, conformational dynamics, and configurational entropy of the antagonist-bound OX1R complex. The current workflow combines standard MD observables such as RMSD, RMSF, radius of gyration, and PCA with configurational-entropy analysis using `PARENT_GPU`.

## Repository Structure

- `orexin1_project.qmd`: main written project report.
- `orexin1_intro_revealjs.qmd`: presentation slides.
- `references.bib`: bibliography used by the Quarto documents.
- `analysis/`: selected derived analysis outputs, summary tables, and publication figures.
- `data/`: raw and prepared structural inputs.
- `scripts/`: helper scripts for trajectory processing, PCA, plotting, and production analysis.
- `system/`: simulation-ready coordinates, topology, and `.mdp` parameter files.
- `visualization/`: PyMOL, VMD, and static visual assets used for structural comparison.

## Included Outputs

This repository intentionally keeps selected summaries and figure-level outputs that support the report, while excluding very large generated trajectory files and other heavy intermediate artifacts from version control.

Included examples:
- equilibrated-system summaries and publication plots
- PCA summaries and rendered figure panels
- production-analysis summaries
- configurational-entropy summary values

Excluded from git:
- production and equilibration run directories under `runs/`
- large trajectory files such as `*.xtc`
- heavy entropy intermediates such as `*.par` and `*.bat`
- Quarto cache and other local generated files

## Reproducibility

The simulation setup is defined in `system/`, and the analysis workflow is captured in `scripts/`. The committed analysis products in `analysis/` are therefore traceable to versioned inputs and scripts even though the heaviest generated files are not stored in the repository.

## Main Documents

- Report: `orexin1_project.qmd`
- Slides: `orexin1_intro_revealjs.qmd`

Both documents can be rendered with Quarto once the required local dependencies are available.
