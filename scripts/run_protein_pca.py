#!/usr/bin/env python3
"""Run backbone PCA on the production trajectory and prepare VMD-friendly outputs."""

from __future__ import annotations

import argparse
import csv
import subprocess
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RUN_DIR = PROJECT_ROOT / "runs" / "production" / "md_results"
DEFAULT_ANALYSIS_DIR = PROJECT_ROOT / "analysis" / "pca"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument("--analysis-dir", type=Path, default=DEFAULT_ANALYSIS_DIR)
    parser.add_argument("--base", default="orexin1_md")
    parser.add_argument("--group", default="Backbone")
    return parser.parse_args()


def run_gmx(cmd: list[str], stdin: str | None = None) -> None:
    subprocess.run(cmd, input=stdin, text=True, check=True, cwd=PROJECT_ROOT)


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


def write_summary(rows: list[dict[str, str | float]], txt_path: Path, csv_path: Path) -> None:
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value", "unit"])
        writer.writeheader()
        writer.writerows(rows)

    lines = ["Backbone PCA summary", ""]
    for row in rows:
        lines.append(f"{row['metric']}: {row['value']} {row['unit']}".rstrip())
    txt_path.write_text("\n".join(lines) + "\n")


def write_vmd_script(path: Path, average_gro: Path, filtered_xtc: Path, combined_xtc: Path) -> None:
    script = f"""mol new {average_gro} type gro waitfor all
mol addfile {filtered_xtc} type xtc waitfor all
mol delrep 0 top
mol representation NewCartoon
mol color Name
mol selection all
mol material AOShiny
mol addrep top

mol new {average_gro} type gro waitfor all
mol addfile {combined_xtc} type xtc waitfor all
mol delrep 0 top
mol representation NewCartoon
mol color ColorID 1
mol selection all
mol material Transparent
mol addrep top
"""
    path.write_text(script)


def main() -> None:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    analysis_dir = args.analysis_dir.resolve()
    analysis_dir.mkdir(parents=True, exist_ok=True)

    base = args.base
    group = args.group

    tpr = run_dir / f"{base}.tpr"
    xtc = run_dir / f"{base}_centered_fin.xtc"
    ndx = run_dir / f"{base}.ndx"

    outputs = {
        "eigenval": analysis_dir / f"{base}_backbone_pca_eigenval.xvg",
        "eigenvec": analysis_dir / f"{base}_backbone_pca_eigenvec.trr",
        "average": analysis_dir / f"{base}_backbone_pca_average.gro",
        "covar_log": analysis_dir / f"{base}_backbone_pca_covar.log",
        "pc1_proj": analysis_dir / f"{base}_backbone_pc1_projection.xvg",
        "pc2_proj": analysis_dir / f"{base}_backbone_pc2_projection.xvg",
        "pc12_2d": analysis_dir / f"{base}_backbone_pc1_pc2_2dproj.xvg",
        "pc1_filtered": analysis_dir / f"{base}_backbone_pc1_filtered.xtc",
        "pc12_filtered": analysis_dir / f"{base}_backbone_pc1_pc2_filtered.xtc",
        "pc1_extremes": analysis_dir / f"{base}_backbone_pc1_extremes.pdb",
        "figure_png": analysis_dir / f"{base}_backbone_pca_grid.png",
        "figure_pdf": analysis_dir / f"{base}_backbone_pca_grid.pdf",
        "figure_svg": analysis_dir / f"{base}_backbone_pca_grid.svg",
        "summary_txt": analysis_dir / f"{base}_backbone_pca_summary.txt",
        "summary_csv": analysis_dir / f"{base}_backbone_pca_summary.csv",
        "vmd_script": analysis_dir / f"{base}_backbone_pca_vmd_load.tcl",
    }

    for path in outputs.values():
        if isinstance(path, Path) and path.exists():
            path.unlink()

    two_groups = f"{group}\n{group}\n"
    three_groups = f"{group}\n{group}\n{group}\n"

    run_gmx(
        [
            "gmx",
            "covar",
            "-s",
            str(tpr),
            "-f",
            str(xtc),
            "-n",
            str(ndx),
            "-o",
            str(outputs["eigenval"]),
            "-v",
            str(outputs["eigenvec"]),
            "-av",
            str(outputs["average"]),
            "-l",
            str(outputs["covar_log"]),
            "-tu",
            "ns",
            "-xvg",
            "none",
        ],
        stdin=two_groups,
    )

    run_gmx(
        [
            "gmx",
            "anaeig",
            "-v",
            str(outputs["eigenvec"]),
            "-eig",
            str(outputs["eigenval"]),
            "-s",
            str(tpr),
            "-f",
            str(xtc),
            "-n",
            str(ndx),
            "-proj",
            str(outputs["pc1_proj"]),
            "-first",
            "1",
            "-last",
            "1",
            "-tu",
            "ns",
            "-xvg",
            "none",
        ],
        stdin=two_groups,
    )
    run_gmx(
        [
            "gmx",
            "anaeig",
            "-v",
            str(outputs["eigenvec"]),
            "-eig",
            str(outputs["eigenval"]),
            "-s",
            str(tpr),
            "-f",
            str(xtc),
            "-n",
            str(ndx),
            "-proj",
            str(outputs["pc2_proj"]),
            "-first",
            "2",
            "-last",
            "2",
            "-tu",
            "ns",
            "-xvg",
            "none",
        ],
        stdin=two_groups,
    )
    run_gmx(
        [
            "gmx",
            "anaeig",
            "-v",
            str(outputs["eigenvec"]),
            "-eig",
            str(outputs["eigenval"]),
            "-s",
            str(tpr),
            "-f",
            str(xtc),
            "-n",
            str(ndx),
            "-2d",
            str(outputs["pc12_2d"]),
            "-first",
            "1",
            "-last",
            "2",
            "-tu",
            "ns",
            "-xvg",
            "none",
        ],
        stdin=two_groups,
    )
    run_gmx(
        [
            "gmx",
            "anaeig",
            "-v",
            str(outputs["eigenvec"]),
            "-eig",
            str(outputs["eigenval"]),
            "-s",
            str(tpr),
            "-f",
            str(xtc),
            "-n",
            str(ndx),
            "-filt",
            str(outputs["pc1_filtered"]),
            "-first",
            "1",
            "-last",
            "1",
            "-tu",
            "ns",
            "-xvg",
            "none",
        ],
        stdin=three_groups,
    )
    run_gmx(
        [
            "gmx",
            "anaeig",
            "-v",
            str(outputs["eigenvec"]),
            "-eig",
            str(outputs["eigenval"]),
            "-s",
            str(tpr),
            "-f",
            str(xtc),
            "-n",
            str(ndx),
            "-filt",
            str(outputs["pc12_filtered"]),
            "-first",
            "1",
            "-last",
            "2",
            "-tu",
            "ns",
            "-xvg",
            "none",
        ],
        stdin=three_groups,
    )
    run_gmx(
        [
            "gmx",
            "anaeig",
            "-v",
            str(outputs["eigenvec"]),
            "-eig",
            str(outputs["eigenval"]),
            "-s",
            str(tpr),
            "-f",
            str(xtc),
            "-n",
            str(ndx),
            "-extr",
            str(outputs["pc1_extremes"]),
            "-first",
            "1",
            "-last",
            "1",
            "-nframes",
            "20",
            "-tu",
            "ns",
            "-xvg",
            "none",
        ],
        stdin=three_groups,
    )

    eigenval = load_xvg(outputs["eigenval"])
    pc1 = load_xvg(outputs["pc1_proj"])
    pc2 = load_xvg(outputs["pc2_proj"])
    pc12 = load_xvg(outputs["pc12_2d"])

    eigenvalues = eigenval[:, 1]
    variance_fraction = eigenvalues / np.sum(eigenvalues)
    cumulative_fraction = np.cumsum(variance_fraction)

    top_n = min(10, len(eigenvalues))
    top_idx = np.arange(1, top_n + 1)

    pc1_time = pc1[:, 0]
    pc1_vals = pc1[:, 1]
    pc2_vals = pc2[:, 1]

    pc1_min_idx = int(np.argmin(pc1_vals))
    pc1_max_idx = int(np.argmax(pc1_vals))
    pc2_min_idx = int(np.argmin(pc2_vals))
    pc2_max_idx = int(np.argmax(pc2_vals))

    summary_rows = [
        {"metric": "Analysis group", "value": group, "unit": ""},
        {"metric": "Frames used", "value": str(len(pc1_time)), "unit": ""},
        {"metric": "PC1 variance explained", "value": f"{100 * variance_fraction[0]:.2f}", "unit": "%"},
        {"metric": "PC2 variance explained", "value": f"{100 * variance_fraction[1]:.2f}", "unit": "%"},
        {"metric": "PC3 variance explained", "value": f"{100 * variance_fraction[2]:.2f}", "unit": "%"},
        {"metric": "PC1+PC2 cumulative variance", "value": f"{100 * cumulative_fraction[1]:.2f}", "unit": "%"},
        {"metric": "PC1+PC2+PC3 cumulative variance", "value": f"{100 * cumulative_fraction[2]:.2f}", "unit": "%"},
        {"metric": "PC1 minimum projection", "value": f"{pc1_vals[pc1_min_idx]:.3f} at {pc1_time[pc1_min_idx]:.1f} ns", "unit": ""},
        {"metric": "PC1 maximum projection", "value": f"{pc1_vals[pc1_max_idx]:.3f} at {pc1_time[pc1_max_idx]:.1f} ns", "unit": ""},
        {"metric": "PC2 minimum projection", "value": f"{pc2_vals[pc2_min_idx]:.3f} at {pc1_time[pc2_min_idx]:.1f} ns", "unit": ""},
        {"metric": "PC2 maximum projection", "value": f"{pc2_vals[pc2_max_idx]:.3f} at {pc1_time[pc2_max_idx]:.1f} ns", "unit": ""},
    ]
    write_summary(summary_rows, outputs["summary_txt"], outputs["summary_csv"])
    write_vmd_script(outputs["vmd_script"], outputs["average"], outputs["pc1_filtered"], outputs["pc12_filtered"])

    mpl.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 600,
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, axes = plt.subplots(2, 2, figsize=(11.8, 9.0))
    axes = axes.ravel()

    axes[0].bar(top_idx, 100 * variance_fraction[:top_n], color="#4C72B0", edgecolor="white", linewidth=0.8)
    axes[0].set_title("Explained Variance By PC")
    axes[0].set_xlabel("Principal component")
    axes[0].set_ylabel("Variance explained (%)")
    axes[0].set_xticks(top_idx)

    cum_n = min(50, len(cumulative_fraction))
    axes[1].plot(np.arange(1, cum_n + 1), 100 * cumulative_fraction[:cum_n], color="#C44E52", lw=2.0)
    axes[1].axhline(80, color="0.4", ls=(0, (3, 2)), lw=1.0)
    axes[1].set_title("Cumulative Variance")
    axes[1].set_xlabel("Principal component")
    axes[1].set_ylabel("Cumulative variance (%)")
    axes[1].set_xlim(1, cum_n)

    axes[2].plot(pc1_time, pc1_vals, color="#55A868", lw=1.7, label=f"PC1 ({100 * variance_fraction[0]:.1f}%)")
    axes[2].plot(pc1_time, pc2_vals, color="#8172B2", lw=1.7, label=f"PC2 ({100 * variance_fraction[1]:.1f}%)")
    axes[2].axhline(0.0, color="0.5", lw=0.9, ls=(0, (2, 2)))
    axes[2].set_title("PC Projections Over Time")
    axes[2].set_xlabel("Time (ns)")
    axes[2].set_ylabel("Projection (nm)")
    axes[2].legend(frameon=False, loc="upper left")

    scatter = axes[3].scatter(pc12[:, 0], pc12[:, 1], c=pc1_time, cmap="viridis", s=18, alpha=0.9, edgecolors="none")
    axes[3].set_title("Conformational Space: PC1 vs PC2")
    axes[3].set_xlabel("PC1 projection (nm)")
    axes[3].set_ylabel("PC2 projection (nm)")
    cbar = fig.colorbar(scatter, ax=axes[3], fraction=0.046, pad=0.04)
    cbar.set_label("Time (ns)")

    for ax in axes:
        ax.grid(True, color="0.93", lw=0.8)

    fig.tight_layout()
    fig.savefig(outputs["figure_png"])
    fig.savefig(outputs["figure_pdf"])
    fig.savefig(outputs["figure_svg"])
    plt.close(fig)


if __name__ == "__main__":
    main()
