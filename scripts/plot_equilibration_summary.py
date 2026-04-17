#!/usr/bin/env python3
"""Create publication-ready plots for the final equilibration stage."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_CSV = PROJECT_ROOT / "analysis" / "6_equilibration_key_observables.csv"
OUTPUT_STEM = PROJECT_ROOT / "analysis" / "6_equilibration_publication_grid"
SUMMARY_TXT = PROJECT_ROOT / "analysis" / "6_equilibration_final_window_summary.txt"
SUMMARY_CSV = PROJECT_ROOT / "analysis" / "6_equilibration_final_window_summary.csv"

FINAL_WINDOW_PS = 100.0
ROLLING_WINDOW_PS = 20.0


PLOT_SPECS = [
    {
        "column": "temperature",
        "title": "Temperature",
        "ylabel": "Temperature (K)",
        "color": "#C44E52",
        "target": 310.15,
    },
    {
        "column": "pressure",
        "title": "Pressure",
        "ylabel": "Pressure (bar)",
        "color": "#4C72B0",
        "target": 1.0,
    },
    {
        "column": "density",
        "title": "Density",
        "ylabel": r"Density (kg m$^{-3}$)",
        "color": "#55A868",
        "target": None,
    },
    {
        "column": "potential",
        "title": "Potential Energy",
        "ylabel": r"Potential (kJ mol$^{-1}$)",
        "color": "#8172B2",
        "target": None,
    },
    {
        "column": "box-x",
        "title": "Box-X",
        "ylabel": "Box-X (nm)",
        "color": "#CCB974",
        "target": None,
    },
    {
        "column": "box-z",
        "title": "Box-Z",
        "ylabel": "Box-Z (nm)",
        "color": "#64B5CD",
        "target": None,
    },
]


def load_csv(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        columns = {name: [] for name in reader.fieldnames or []}
        for row in reader:
            for name, value in row.items():
                columns[name].append(float(value))
    return {name: np.asarray(values, dtype=float) for name, values in columns.items()}


def rolling_mean(values: np.ndarray, window_points: int) -> np.ndarray:
    if window_points <= 1:
        return values.copy()
    if window_points > len(values):
        window_points = len(values)
    kernel = np.ones(window_points, dtype=float) / window_points
    padded = np.pad(values, (window_points // 2, window_points - 1 - window_points // 2), mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def format_stat(mean: float, std: float, unit: str) -> str:
    if abs(mean) >= 1.0e4:
        return f"{mean:,.1f} ± {std:,.1f} {unit}".replace(",", "")
    if abs(mean) >= 100:
        return f"{mean:.2f} ± {std:.2f} {unit}"
    return f"{mean:.3f} ± {std:.3f} {unit}"


def main() -> None:
    data = load_csv(INPUT_CSV)
    time_ps = data["time_ps"]

    dt_ps = np.median(np.diff(time_ps))
    rolling_points = max(1, int(round(ROLLING_WINDOW_PS / dt_ps)))
    final_mask = time_ps >= (time_ps.max() - FINAL_WINDOW_PS)

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
            "legend.fontsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, axes = plt.subplots(2, 3, figsize=(15.0, 8.8))
    axes = axes.ravel()

    summary_rows: list[dict[str, str | float]] = []

    for ax, spec in zip(axes, PLOT_SPECS):
        values = data[spec["column"]]
        smooth = rolling_mean(values, rolling_points)

        raw_line = ax.plot(time_ps, values, color="0.80", lw=0.9, alpha=0.95, label="Raw")[0]
        smooth_line = ax.plot(
            time_ps,
            smooth,
            color=spec["color"],
            lw=2.0,
            label=f"{int(ROLLING_WINDOW_PS)} ps rolling mean",
        )[0]

        ax.axvspan(time_ps.max() - FINAL_WINDOW_PS, time_ps.max(), color="#EAEAF2", alpha=0.9, zorder=0)

        if spec["target"] is not None:
            ax.axhline(spec["target"], color="0.25", lw=1.0, ls=(0, (4, 2)), label="Target")

        final_values = values[final_mask]
        final_mean = float(np.mean(final_values))
        final_std = float(np.std(final_values, ddof=1))
        last_value = float(values[-1])

        ax.axhline(final_mean, color=spec["color"], lw=1.2, ls=(0, (1.5, 2.0)), alpha=0.95, label="Final 100 ps mean")

        annotation = f"Final 100 ps: {final_mean:.2f} ± {final_std:.2f}\nLast frame: {last_value:.2f}"
        ax.text(
            0.02,
            0.96,
            annotation,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.5,
            bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "0.85", "alpha": 0.95},
        )

        ax.set_title(spec["title"], pad=8)
        ax.set_ylabel(spec["ylabel"])
        ax.set_xlabel("Time (ps)")
        ax.grid(True, color="0.92", lw=0.8)
        ax.set_xlim(time_ps.min(), time_ps.max())

        summary_rows.append(
            {
                "observable": spec["title"],
                "column": spec["column"],
                "final_window_ps": FINAL_WINDOW_PS,
                "final_window_mean": final_mean,
                "final_window_std": final_std,
                "last_frame_value": last_value,
                "target": "" if spec["target"] is None else spec["target"],
            }
        )

    fig.tight_layout()

    for suffix in (".png", ".pdf", ".svg"):
        fig.savefig(OUTPUT_STEM.with_suffix(suffix))
    plt.close(fig)

    with SUMMARY_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)

    lines = [
        "Final equilibration window summary",
        f"Source: {INPUT_CSV}",
        f"Window: last {FINAL_WINDOW_PS:.0f} ps of the 500 ps stage-6 equilibration",
        "",
    ]
    units = {
        "Temperature": "K",
        "Pressure": "bar",
        "Density": "kg/m^3",
        "Potential Energy": "kJ/mol",
        "Box-X": "nm",
        "Box-Z": "nm",
    }
    for row in summary_rows:
        title = str(row["observable"])
        unit = units[title]
        target = row["target"]
        lines.append(f"{title}:")
        lines.append(
            f"  final 100 ps mean ± SD = {format_stat(float(row['final_window_mean']), float(row['final_window_std']), unit)}"
        )
        lines.append(f"  last frame = {float(row['last_frame_value']):.3f} {unit}")
        if target != "":
            lines.append(f"  target = {float(target):.3f} {unit}")
        lines.append("")

    SUMMARY_TXT.write_text("\n".join(lines))


if __name__ == "__main__":
    main()
