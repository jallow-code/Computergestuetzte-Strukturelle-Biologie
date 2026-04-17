# Computergestuetzte Strukturelle Biologie

This repository contains a molecular dynamics project on the human orexin-1 receptor (OX1R) bound to the selective antagonist GSK1059865. It combines simulation setup, selected analysis outputs, visualization assets, and Quarto-based written material for the final report and slides.

## Repository Structure

- `reports/`: Quarto source documents for the main project report and course report.
- `slides/`: Quarto reveal.js presentation source.
- `references/`: shared bibliography files.
- `styles/`: CSL and document styling assets used by the reports.
- `analysis/`: selected figures, summaries, and reduced analysis outputs tracked in git.
- `data/`: raw structures and prepared ligand inputs.
- `scripts/`: analysis and workflow helper scripts.
- `system/`: coordinates, topology, and GROMACS parameter files for the MD system.
- `visualization/`: PyMOL, VMD, and image assets for structural inspection.

## Main Documents

- Main report: `reports/orexin1_project.qmd`
- Coursework report: `reports/Exercises_report.qmd`
- Slides: `slides/orexin1_intro_revealjs.qmd`

## Scope of Versioned Data

This repository tracks the project inputs, scripts, selected reduced outputs, and publication-ready figures. Large generated artifacts such as raw trajectories, cached render outputs, intermediate entropy files, and temporary local files are intentionally excluded from version control.

## Reproducibility

The MD system definition lives in `system/`, the prepared structural inputs live in `data/`, and the workflow logic lives in `scripts/`. The tracked files in `analysis/` are therefore tied to versioned inputs and analysis code even though the heaviest run products under `runs/` are not committed.
