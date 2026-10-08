#!/usr/bin/env python3
"""Compare regenerated 2609 phonon bands with the immutable retry reference."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml


def load_band(path: Path) -> dict:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    segments = []
    offset = 0
    for count in data["segment_nqpoint"]:
        segment = data["phonon"][offset:offset + count]
        offset += count
        segments.append({
            "distance": np.asarray([point["distance"] for point in segment], dtype=float),
            "frequency": np.asarray([
                [band["frequency"] for band in point["band"]] for point in segment
            ], dtype=float),
        })
    return {"raw": data, "segments": segments}


def compare_system(reference_path: Path, current_path: Path) -> dict:
    reference = load_band(reference_path)
    current = load_band(current_path)
    if len(reference["segments"]) != len(current["segments"]):
        raise ValueError("Band-path segment counts differ")

    errors = []
    ref_values = []
    current_values = []
    segment_rows = []
    for index, (ref_segment, current_segment) in enumerate(zip(reference["segments"], current["segments"]), 1):
        if ref_segment["frequency"].shape != current_segment["frequency"].shape:
            raise ValueError(
                f"Segment {index} band shapes differ: {ref_segment['frequency'].shape} vs "
                f"{current_segment['frequency'].shape}"
            )
        delta = current_segment["frequency"] - ref_segment["frequency"]
        errors.append(delta.ravel())
        ref_values.append(ref_segment["frequency"].ravel())
        current_values.append(current_segment["frequency"].ravel())
        segment_rows.append({
            "segment": index,
            "rmse_thz": float(np.sqrt(np.mean(delta ** 2))),
            "mae_thz": float(np.mean(np.abs(delta))),
            "max_abs_thz": float(np.max(np.abs(delta))),
        })
    error = np.concatenate(errors)
    reference_frequency = np.concatenate(ref_values)
    current_frequency = np.concatenate(current_values)
    correlation = float(np.corrcoef(reference_frequency, current_frequency)[0, 1])

    window_counts = {}
    for lower, upper in ((0, 12), (35, 46), (60, 70)):
        name = f"{lower}-{upper}_THz"
        window_counts[name] = {
            "retry": int(np.count_nonzero((reference_frequency >= lower) & (reference_frequency <= upper))),
            "2609": int(np.count_nonzero((current_frequency >= lower) & (current_frequency <= upper))),
        }
    return {
        "reference": str(reference_path),
        "current": str(current_path),
        "frequency_correlation": correlation,
        "rmse_thz": float(np.sqrt(np.mean(error ** 2))),
        "mae_thz": float(np.mean(np.abs(error))),
        "max_abs_thz": float(np.max(np.abs(error))),
        "segment_metrics": segment_rows,
        "frequency_window_point_counts": window_counts,
        "reference_data": reference,
        "current_data": current,
    }


def plot_comparison(system: str, result: dict, output: Path, dpi: int) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), sharey=True)
    for axis, key, title, color in (
        (axes[0], "reference_data", "phonon_retry reference", "#666666"),
        (axes[1], "current_data", "phonon_retry_2609", "#8C2D04"),
    ):
        for segment in result[key]["segments"]:
            distance = segment["distance"]
            for branch in segment["frequency"].T:
                axis.plot(distance, branch, color=color, linewidth=0.45, alpha=0.75)
        axis.axhspan(35, 46, color="#24B486", alpha=0.08)
        axis.axhspan(60, 70, color="#8D76E2", alpha=0.08)
        axis.set_title(title)
        axis.set_xlabel("Wave-vector path")
        axis.grid(alpha=0.2)
    axes[0].set_ylabel("Frequency (THz)")
    fig.suptitle(
        f"{system}: retry vs 2609 · RMSE {result['rmse_thz']:.4f} THz · "
        f"r = {result['frequency_correlation']:.6f}"
    )
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, default=Path("reference_retry/output"))
    parser.add_argument("--current-dir", type=Path, default=Path("output_2609"))
    parser.add_argument("--outdir", type=Path, default=Path("comparison"))
    parser.add_argument("--systems", nargs="+", default=["CA", "LA", "NEt4-CA"])
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    report = {"systems": {}}
    for system in args.systems:
        result = compare_system(
            args.reference_dir / f"{system}_band.yaml",
            args.current_dir / f"{system}_band.yaml",
        )
        plot_comparison(system, result, args.outdir / f"{system}_retry_vs_2609.png", args.dpi)
        report["systems"][system] = {
            key: value for key, value in result.items()
            if key not in {"reference_data", "current_data"}
        }
        print(
            f"{system}: RMSE={result['rmse_thz']:.6f} THz, "
            f"MAE={result['mae_thz']:.6f} THz, r={result['frequency_correlation']:.8f}"
        )

    report_path = args.outdir / "retry_vs_2609_metrics.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Saved {report_path}")


if __name__ == "__main__":
    main()
