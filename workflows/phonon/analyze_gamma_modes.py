#!/usr/bin/env python3
"""Rank persisted Gamma modes by framework, azide, and organic participation."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


DEFAULT_WINDOWS = {
    "low_frequency": (0.1, 12.0),
    "organic_intermediate": (35.0, 46.0),
    "azide_stretch": (60.0, 70.0),
}


def load_modes(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as archive:
        data = {key: archive[key] for key in archive.files}
    group_names = [str(value) for value in data["group_names"]]
    participation = data.get("mass_weighted_participation", data["participation"])
    data["participation_by_group"] = {
        name: participation[:, index]
        for index, name in enumerate(group_names)
    }
    return data


def ranked_modes(data: dict, group: str, lower: float, upper: float, count: int) -> list[dict]:
    frequencies = data["frequencies_thz"]
    values = data["participation_by_group"][group]
    candidates = np.flatnonzero((frequencies >= lower) & (frequencies <= upper))
    order = candidates[np.argsort(values[candidates])[::-1]][:count]
    rows = []
    for mode_index in order:
        rows.append({
            "mode_index_1based": int(mode_index + 1),
            "frequency_thz": float(frequencies[mode_index]),
            "target_group": group,
            "target_participation": float(values[mode_index]),
            "imaginary_fraction": float(data["imaginary_fraction"][mode_index]),
            "participation": {
                name: float(group_values[mode_index])
                for name, group_values in data["participation_by_group"].items()
            },
        })
    return rows


def write_csv(path: Path, systems: list[dict]) -> None:
    group_names = sorted({
        group
        for system in systems
        for selection in system["selections"].values()
        for row in selection
        for group in row["participation"]
    })
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "system", "window", "mode_index_1based", "frequency_thz",
            "target_group", "target_participation", "imaginary_fraction",
            *[f"{group}_participation" for group in group_names],
        ])
        for system in systems:
            for window, rows in system["selections"].items():
                for row in rows:
                    writer.writerow([
                        system["system"], window, row["mode_index_1based"],
                        f"{row['frequency_thz']:.8f}", row["target_group"],
                        f"{row['target_participation']:.8f}",
                        f"{row['imaginary_fraction']:.3e}",
                        *[f"{row['participation'].get(group, 0.0):.8f}" for group in group_names],
                    ])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archives", type=Path, nargs="+", help="*_gamma_modes.npz files")
    parser.add_argument("--top", type=int, default=8, help="Candidates per window")
    parser.add_argument("--output", type=Path, default=Path("output/gamma_mode_candidates.json"))
    args = parser.parse_args()

    systems = []
    for archive in args.archives:
        data = load_modes(archive)
        system = archive.name.removesuffix("_gamma_modes.npz")
        selections = {
            "azide_stretch": ranked_modes(
                data, "azide", *DEFAULT_WINDOWS["azide_stretch"], args.top
            )
        }
        if "organic_cation" in data["participation_by_group"]:
            selections = {
                "low_frequency_organic": ranked_modes(
                    data, "organic_cation", *DEFAULT_WINDOWS["low_frequency"], args.top
                ),
                "intermediate_organic": ranked_modes(
                    data, "organic_cation", *DEFAULT_WINDOWS["organic_intermediate"], args.top
                ),
                **selections,
            }
        systems.append({"system": system, "archive": str(archive), "selections": selections})

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({
        "note": (
            "Participation ranks candidate modes but does not establish torsion, libration, "
            "deformation, or stretching character; inspect rendered eigenvectors before labeling."
        ),
        "windows_thz": {name: list(bounds) for name, bounds in DEFAULT_WINDOWS.items()},
        "systems": systems,
    }, indent=2) + "\n", encoding="utf-8")
    csv_path = args.output.with_suffix(".csv")
    write_csv(csv_path, systems)
    print(f"Saved {args.output}")
    print(f"Saved {csv_path}")


if __name__ == "__main__":
    main()
