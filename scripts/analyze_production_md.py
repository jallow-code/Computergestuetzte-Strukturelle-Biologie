#!/usr/bin/env python3
"""Run standard production MD analysis for the protein-ligand complex."""

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
DEFAULT_ANALYSIS_DIR = PROJECT_ROOT / "analysis" / "production"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument("--base", default="orexin1_md")
    parser.add_argument("--analysis-dir", type=Path, default=DEFAULT_ANALYSIS_DIR)
    return parser.parse_args()


def run_gmx(cmd: list[str], stdin: str | None = None) -> None:
    subprocess.run(
        cmd,
        input=stdin,
        text=True,
        check=True,
        cwd=PROJECT_ROOT,
    )


def load_xvg(path: Path) -> np.ndarray:
    rows: list[list[float]] = []
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", "@")):
                continue
            rows.append([float(token) for token in stripped.split()])
    return np.asarray(rows, dtype=float)


def rolling_mean(values: np.ndarray, window_points: int) -> np.ndarray:
    if window_points <= 1:
        return values.copy()
    if window_points > len(values):
        window_points = len(values)
    kernel = np.ones(window_points, dtype=float) / window_points
    padded = np.pad(values, (window_points // 2, window_points - 1 - window_points // 2), mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def summarize_time_series(name: str, values: np.ndarray, unit: str) -> dict[str, str | float]:
    return {
        "metric": name,
        "mean": float(np.mean(values)),
        "std": float(np.std(values, ddof=1)),
        "min": float(np.min(values)),
        "max": float(np.max(values)),
        "last": float(values[-1]),
        "unit": unit,
    }


def write_summary(summary_rows: list[dict[str, str | float]], path_txt: Path, path_csv: Path) -> None:
    with path_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "mean", "std", "min", "max", "last", "unit"])
        writer.writeheader()
        writer.writerows(summary_rows)

    lines = ["Production MD summary", ""]
    for row in summary_rows:
        lines.append(f"{row['metric']}:")
        lines.append(
            f"  mean ± SD = {float(row['mean']):.3f} ± {float(row['std']):.3f} {row['unit']}"
        )
        lines.append(f"  range = {float(row['min']):.3f} to {float(row['max']):.3f} {row['unit']}")
        lines.append(f"  last frame = {float(row['last']):.3f} {row['unit']}")
        lines.append("")
    path_txt.write_text("\n".join(lines))


def main() -> None:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    analysis_dir = args.analysis_dir.resolve()
    analysis_dir.mkdir(parents=True, exist_ok=True)

    base = args.base
    tpr = run_dir / f"{base}.tpr"
    ndx = run_dir / f"{base}.ndx"
    xtc = run_dir / f"{base}_centered_fin.xtc"

    outputs = {
        "protein_rmsd": analysis_dir / f"{base}_protein_backbone_rmsd.xvg",
        "ligand_rmsd": analysis_dir / f"{base}_ligand_rmsd_vs_proteinfit.xvg",
        "calpha_rmsf": analysis_dir / f"{base}_calpha_rmsf.xvg",
        "rg": analysis_dir / f"{base}_protein_rg.xvg",
        "com_distance": analysis_dir / f"{base}_protein_ligand_com_distance.xvg",
        "hbond": analysis_dir / f"{base}_protein_ligand_hbonds.xvg",
        "hbond_index": analysis_dir / f"{base}_protein_ligand_hbonds.ndx",
        "min_distance": analysis_dir / f"{base}_protein_ligand_min_distance.xvg",
        "summary_txt": analysis_dir / f"{base}_analysis_summary.txt",
        "summary_csv": analysis_dir / f"{base}_analysis_summary.csv",
        "figure_png": analysis_dir / f"{base}_analysis_grid.png",
        "figure_pdf": analysis_dir / f"{base}_analysis_grid.pdf",
        "figure_svg": analysis_dir / f"{base}_analysis_grid.svg",
    }

    for output_path in outputs.values():
        if isinstance(output_path, Path) and output_path.exists():
            output_path.unlink()

    run_gmx(
        [
            "gmx",
            "rms",
            "-s",
            str(tpr),
            "-f",
            str(xtc),
            "-n",
            str(ndx),
            "-o",
            str(outputs["protein_rmsd"]),
            "-tu",
            "ns",
            "-xvg",
            "none",
        ],
        stdin="Backbone\nBackbone\n",
    )
    run_gmx(
        [
            "gmx",
            "rms",
            "-s",
            str(tpr),
            "-f",
            str(xtc),
            "-n",
            str(ndx),
            "-o",
            str(outputs["ligand_rmsd"]),
            "-tu",
            "ns",
            "-xvg",
            "none",
        ],
        stdin="Backbone\nNRE\n",
    )
    run_gmx(
        [
            "gmx",
            "rmsf",
            "-s",
            str(tpr),
            "-f",
            str(xtc),
            "-n",
            str(ndx),
            "-o",
            str(outputs["calpha_rmsf"]),
            "-res",
            "-xvg",
            "none",
        ],
        stdin="C-alpha\n",
    )
    run_gmx(
        [
            "gmx",
            "gyrate",
            "-s",
            str(tpr),
            "-f",
            str(xtc),
            "-n",
            str(ndx),
            "-o",
            str(outputs["rg"]),
            "-tu",
            "ns",
            "-xvg",
            "none",
        ],
        stdin="Protein\n",
    )
    run_gmx(
        [
            "gmx",
            "distance",
            "-s",
            str(tpr),
            "-f",
            str(xtc),
            "-n",
            str(ndx),
            "-select",
            'com of group "Protein" plus com of group "NRE"',
            "-oall",
            str(outputs["com_distance"]),
            "-tu",
            "ns",
            "-xvg",
            "none",
        ]
    )
    run_gmx(
        [
            "gmx",
            "hbond",
            "-s",
            str(tpr),
            "-f",
            str(xtc),
            "-n",
            str(ndx),
            "-r",
            'group "Protein"',
            "-t",
            'group "NRE"',
            "-o",
            str(outputs["hbond_index"]),
            "-num",
            str(outputs["hbond"]),
            "-tu",
            "ns",
            "-xvg",
            "none",
        ]
    )
    run_gmx(
        [
            "gmx",
            "pairdist",
            "-s",
            str(tpr),
            "-f",
            str(xtc),
            "-n",
            str(ndx),
            "-ref",
            'group "NRE"',
            "-sel",
            'group "Protein"',
            "-o",
            str(outputs["min_distance"]),
            "-tu",
            "ns",
            "-xvg",
            "none",
        ]
    )

    protein_rmsd = load_xvg(outputs["protein_rmsd"])
    ligand_rmsd = load_xvg(outputs["ligand_rmsd"])
    calpha_rmsf = load_xvg(outputs["calpha_rmsf"])
    rg = load_xvg(outputs["rg"])
    com_distance = load_xvg(outputs["com_distance"])
    hbonds = load_xvg(outputs["hbond"])
    min_distance = load_xvg(outputs["min_distance"])

    time_ns = protein_rmsd[:, 0]
    dt_ns = np.median(np.diff(time_ns))
    smooth_points = max(1, int(round(2.0 / dt_ns)))

    summary_rows = [
        summarize_time_series("Protein backbone RMSD", protein_rmsd[:, 1], "nm"),
        summarize_time_series("Ligand RMSD (protein-fit)", ligand_rmsd[:, 1], "nm"),
        summarize_time_series("Protein radius of gyration", rg[:, 1], "nm"),
        summarize_time_series("C-alpha RMSF", calpha_rmsf[:, 1], "nm"),
    ]
    write_summary(summary_rows, outputs["summary_txt"], outputs["summary_csv"])

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
            "legend.fontsize": 8.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.6))
    axes = axes.ravel()

    protein_avg = float(np.mean(protein_rmsd[:, 1]))
    ligand_avg = float(np.mean(ligand_rmsd[:, 1]))
    axes[0].plot(time_ns, protein_rmsd[:, 1], color="#C44E52", lw=1.8, label=f"Protein backbone ({protein_avg:.3f} nm)")
    axes[0].plot(time_ns, ligand_rmsd[:, 1], color="#4C72B0", lw=1.8, label=f"Ligand ({ligand_avg:.3f} nm)")
    axes[0].plot(time_ns, rolling_mean(protein_rmsd[:, 1], smooth_points), color="#C44E52", lw=0.9, alpha=0.35)
    axes[0].plot(time_ns, rolling_mean(ligand_rmsd[:, 1], smooth_points), color="#4C72B0", lw=0.9, alpha=0.35)
    axes[0].set_title("RMSD")
    axes[0].set_xlabel("Time (ns)")
    axes[0].set_ylabel("RMSD (nm)")
    axes[0].legend(loc="upper left", frameon=True, framealpha=0.9)

    axes[1].plot(calpha_rmsf[:, 0], calpha_rmsf[:, 1], color="#55A868", lw=1.7)
    axes[1].fill_between(calpha_rmsf[:, 0], calpha_rmsf[:, 1], 0, color="#55A868", alpha=0.18)
    axes[1].set_title("C-alpha RMSF")
    axes[1].set_xlabel("Residue")
    axes[1].set_ylabel("RMSF (nm)")

    axes[2].plot(rg[:, 0], rg[:, 1], color="#8172B2", lw=1.8)
    axes[2].axhline(np.mean(rg[:, 1]), color="#8172B2", ls=(0, (2, 2)), lw=1.1, alpha=0.8)
    axes[2].set_title("Protein Radius Of Gyration")
    axes[2].set_xlabel("Time (ns)")
    axes[2].set_ylabel("Rg (nm)")

    for ax in axes:
        ax.grid(True, color="0.92", lw=0.8)
        ax.set_xlim(left=0)

    fig.tight_layout()
    fig.savefig(outputs["figure_png"])
    fig.savefig(outputs["figure_pdf"])
    fig.savefig(outputs["figure_svg"])
    plt.close(fig)


if __name__ == "__main__":
    main()
