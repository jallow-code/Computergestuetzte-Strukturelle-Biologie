#!/usr/bin/env python3
"""Generate publication-quality plots for MD analysis."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = PROJECT_ROOT / "analysis" / "production"
PCA_DIR = PROJECT_ROOT / "analysis" / "pca"
OUTPUT_DIR = PROJECT_ROOT / "analysis" / "publication_figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_xvg(path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", "@")):
                continue
            if stripped == "&":
                continue
            rows.append([float(token) for token in stripped.split()])
    return np.asarray(rows, dtype=float)


def set_pub_style():
    """Set aesthetic parameters for publication-quality plots."""
    sns.set_theme(style="ticks", context="paper")
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "figure.dpi": 300,
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def plot_rmsd():
    protein_rmsd = load_xvg(ANALYSIS_DIR / "orexin1_md_protein_backbone_rmsd.xvg")
    ligand_rmsd = load_xvg(ANALYSIS_DIR / "orexin1_md_ligand_rmsd_vs_proteinfit.xvg")

    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.plot(protein_rmsd[:, 0], protein_rmsd[:, 1], label="Protein Backbone", color="#4C72B0", alpha=0.8, lw=1.2)
    ax.plot(ligand_rmsd[:, 0], ligand_rmsd[:, 1], label="Ligand (GSK1059865)", color="#C44E52", alpha=0.8, lw=1.2)

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("RMSD (nm)")
    ax.set_xlim(0, 100)
    ax.legend(frameon=False)
    sns.despine()
    
    fig.savefig(OUTPUT_DIR / "figure_rmsd.pdf")
    fig.savefig(OUTPUT_DIR / "figure_rmsd.png")
    plt.close(fig)


def plot_rmsf():
    rmsf = load_xvg(ANALYSIS_DIR / "orexin1_md_calpha_rmsf.xvg")

    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.fill_between(rmsf[:, 0], rmsf[:, 1], 0, color="#55A868", alpha=0.2)
    ax.plot(rmsf[:, 0], rmsf[:, 1], color="#55A868", lw=1.0)

    ax.set_xlabel("Residue Number")
    ax.set_ylabel("RMSF (nm)")
    ax.set_xlim(rmsf[0, 0], rmsf[-1, 0])
    sns.despine()
    
    fig.savefig(OUTPUT_DIR / "figure_rmsf.pdf")
    fig.savefig(OUTPUT_DIR / "figure_rmsf.png")
    plt.close(fig)


def plot_pca():
    pca_2d = load_xvg(PCA_DIR / "orexin1_md_backbone_pc1_pc2_2dproj.xvg")
    time = np.linspace(0, 100, len(pca_2d))

    fig, ax = plt.subplots(figsize=(5, 4))
    scatter = ax.scatter(pca_2d[:, 0], pca_2d[:, 1], c=time, cmap="viridis", s=10, alpha=0.6, edgecolors="none")
    
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("Time (ns)")
    
    ax.set_xlabel("PC1 (nm)")
    ax.set_ylabel("PC2 (nm)")
    sns.despine()
    
    fig.savefig(OUTPUT_DIR / "figure_pca_2d.pdf")
    fig.savefig(OUTPUT_DIR / "figure_pca_2d.png")
    plt.close(fig)


def plot_activation_markers():
    distances = load_xvg(ANALYSIS_DIR / "activation_distances.xvg")
    
    fig, ax = plt.subplots(figsize=(5, 3.5))
    ax.plot(distances[:, 0], distances[:, 1], label="3.50 - 6.30", color="#4C72B0", lw=1.2)
    ax.plot(distances[:, 0], distances[:, 2], label="3.50 - 7.53", color="#C44E52", lw=1.2)
    ax.plot(distances[:, 0], distances[:, 3], label="6.30 - 7.53", color="#55A868", lw=1.2)

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Distance (nm)")
    ax.set_xlim(0, 100)
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.15), ncol=3)
    sns.despine()
    
    fig.savefig(OUTPUT_DIR / "figure_activation_markers.pdf")
    fig.savefig(OUTPUT_DIR / "figure_activation_markers.png")
    plt.close(fig)


def main():
    set_pub_style()
    plot_rmsd()
    plot_rmsf()
    plot_pca()
    plot_activation_markers()
    print(f"Publication-ready plots generated in {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
