#!/usr/bin/env python3
"""Render persisted Phonopy Gamma modes as consistent MatterVis PNG/GIF media."""

from __future__ import annotations

import argparse
import io
import json
import os
import sys
from pathlib import Path

import gemmi
import numpy as np
from ase import Atoms
from ase.io import write as ase_write
from ase.neighborlist import neighbor_list
from PIL import Image, ImageDraw, ImageFont

REFERENCE_ROOT = Path(os.environ.get(
    "CA_SERIES_MATTERVIS_ROOT",
    "",
))
MATTERVIS = REFERENCE_ROOT / "algorithm" / "MatterVis" / "algorithm" / "MatterVis-world-vectors"
if not MATTERVIS.is_dir():
    raise ImportError(
        "MatterVis is an optional external renderer; set CA_SERIES_MATTERVIS_ROOT "
        "to its local checkout before using render_gamma_modes.py"
    )
if str(MATTERVIS) not in sys.path:
    sys.path.insert(0, str(MATTERVIS))

from crystal_viewer.renderer import build_figure
from crystal_viewer.scene import build_scene_from_atoms, scene_ops, scene_style


ELEMENT_COLORS = {
    "Cu": "#B87333",
    "Pb": "#5B6F8E",
    "N": "#377EB8",
    "C": "#555555",
    "H": "#D9D9D9",
}
ARROW_COLOR = "crimson"


def camera_dict(eye: tuple[float, float, float]) -> dict:
    return {
        "eye": {"x": eye[0], "y": eye[1], "z": eye[2]},
        "center": {"x": 0.0, "y": 0.0, "z": 0.0},
        "up": {"x": 0.0, "y": 0.0, "z": 1.0},
        "projection": {"type": "orthographic"},
    }


def axis_camera(matrix: np.ndarray, axis: str, eye_distance: float = 1.8) -> dict:
    real = [np.asarray(matrix[index], dtype=float) for index in range(3)]
    reciprocal_matrix = np.linalg.inv(np.asarray(matrix, dtype=float))
    reciprocal = [reciprocal_matrix[:, index] for index in range(3)]
    vectors = {
        key: vector / np.linalg.norm(vector)
        for key, vector in zip(("a", "b", "c", "a*", "b*", "c*"), real + reciprocal)
    }
    view = vectors[axis]
    up_key = {"a": "c", "b": "c", "c": "b", "a*": "c*", "b*": "c*", "c*": "b*"}[axis]
    up = vectors[up_key] - np.dot(vectors[up_key], view) * view
    up /= np.linalg.norm(up)
    return camera_dict(tuple(eye_distance * view)) | {
        "up": {"x": float(up[0]), "y": float(up[1]), "z": float(up[2])}
    }


def load_archive(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key] for key in archive.files}


def build_atoms(data: dict) -> Atoms:
    return Atoms(
        symbols=[str(symbol) for symbol in data["symbols"]],
        cell=data["cell"],
        scaled_positions=data["scaled_positions"],
        pbc=True,
    )


def atom_groups(data: dict) -> tuple[list[str], dict[str, np.ndarray]]:
    names = [str(value) for value in data["group_names"]]
    values = data.get("mass_weighted_participation", data["participation"])
    participation = {
        name: values[:, index] for index, name in enumerate(names)
    }
    return names, participation


def payload_from_atoms(atoms: Atoms) -> tuple[list[dict], gemmi.UnitCell, np.ndarray]:
    matrix = np.asarray(atoms.cell, dtype=float)
    lengths = atoms.cell.lengths()
    angles = atoms.cell.angles()
    cell = gemmi.UnitCell(*lengths, *angles)
    payload = []
    for index, (symbol, position, fractional) in enumerate(zip(
        atoms.get_chemical_symbols(), atoms.get_positions(), atoms.get_scaled_positions()
    )):
        payload.append({
            "label": f"{symbol}{index + 1}",
            "elem": symbol,
            "cart": np.asarray(position, dtype=float),
            "frac": np.asarray(fractional, dtype=float),
            "occ": 1.0,
            "dg": ".",
            "da": ".",
            "_source_index": index,
        })
    return payload, cell, matrix


def classified_source_indices(atoms: Atoms) -> dict[str, np.ndarray]:
    """Return source indices for organic cations, azides, and framework."""
    symbols = np.asarray(atoms.get_chemical_symbols())
    left, right, distances = neighbor_list("ijd", atoms, 1.9)
    organic_n = {
        int(i if symbols[i] == "N" else j)
        for i, j, distance in zip(left, right, distances)
        if {symbols[i], symbols[j]} == {"C", "N"} and distance < 1.75
    }
    organic = np.asarray([
        index for index, symbol in enumerate(symbols)
        if symbol in {"C", "H"} or index in organic_n
    ], dtype=int)
    azide = np.asarray([
        index for index, symbol in enumerate(symbols)
        if symbol == "N" and index not in organic_n
    ], dtype=int)
    assigned = set(organic) | set(azide)
    framework = np.asarray([
        index for index in range(len(atoms)) if index not in assigned
    ], dtype=int)
    return {"all": np.arange(len(atoms)), "organic": organic, "azide": azide, "framework": framework}


def focused_payload(
    atoms: Atoms,
    source_indices: np.ndarray,
) -> tuple[list[dict], gemmi.UnitCell, np.ndarray, list[tuple[int, int]]]:
    """Build an unwrapped, molecule-readable cluster for selected source atoms."""
    symbols = np.asarray(atoms.get_chemical_symbols())
    positions = np.asarray(atoms.get_positions(), dtype=float)
    selected = set(int(index) for index in source_indices)
    left, right, distances = neighbor_list("ijd", atoms, 3.2)

    adjacency = {index: [] for index in selected}
    bond_pairs_source = set()
    for i, j, distance in zip(left, right, distances):
        i, j = int(i), int(j)
        if i not in selected or j not in selected or i == j:
            continue
        pair = {symbols[i], symbols[j]}
        bonded = (
            (pair == {"C", "H"} and distance < 1.25)
            or (pair == {"C", "N"} and distance < 1.75)
            or (pair == {"C"} and distance < 1.80)
            or (pair == {"N"} and distance < 1.55)
        )
        if bonded:
            adjacency[i].append(j)
            bond_pairs_source.add(tuple(sorted((i, j))))

    components = []
    remaining = set(selected)
    while remaining:
        root = remaining.pop()
        component = {root}
        stack = [root]
        while stack:
            for neighbour in adjacency[stack.pop()]:
                if neighbour in remaining:
                    remaining.remove(neighbour)
                    component.add(neighbour)
                    stack.append(neighbour)
        components.append(sorted(component))

    unwrapped = {}
    for component in components:
        root = component[0]
        unwrapped[root] = positions[root]
        stack = [root]
        visited = {root}
        while stack:
            current = stack.pop()
            for neighbour in adjacency[current]:
                if neighbour in visited:
                    continue
                vector = atoms.get_distance(current, neighbour, mic=True, vector=True)
                unwrapped[neighbour] = unwrapped[current] + vector
                visited.add(neighbour)
                stack.append(neighbour)
        for index in component:
            unwrapped.setdefault(index, positions[index])

    component_centers = [np.mean([unwrapped[index] for index in component], axis=0) for component in components]
    global_center = np.mean(component_centers, axis=0)
    centered_positions = {
        index: unwrapped[index] - global_center for component in components for index in component
    }
    extents = np.ptp(np.vstack(list(centered_positions.values())), axis=0)
    box = max(18.0, float(extents.max()) + 6.0)
    matrix = np.eye(3) * box
    cell = gemmi.UnitCell(box, box, box, 90.0, 90.0, 90.0)

    ordered = sorted(selected)
    payload = []
    for source in ordered:
        position = centered_positions[source]
        payload.append({
            "label": f"{symbols[source]}{source + 1}",
            "elem": str(symbols[source]),
            "cart": position,
            "frac": position / box + 0.5,
            "occ": 1.0,
            "dg": ".",
            "da": ".",
            "_source_index": source,
        })
    return payload, cell, matrix, sorted(bond_pairs_source)


def mode_arrows(
    positions: np.ndarray,
    symbols: list[str],
    mode: np.ndarray,
    *,
    length: float,
    minimum_relative_amplitude: float,
    heavy_atoms_only: bool,
    source_indices: np.ndarray | None = None,
) -> tuple[list[dict], dict]:
    amplitudes = np.linalg.norm(mode, axis=1)
    maximum = float(amplitudes.max())
    relative = amplitudes / maximum
    selected = [
        index for index, symbol in enumerate(symbols)
        if relative[index] >= minimum_relative_amplitude
        and (not heavy_atoms_only or symbol != "H")
    ]
    arrows = []
    if source_indices is None:
        source_indices = np.arange(len(positions), dtype=int)
    for index in selected:
        source_index = int(source_indices[index])
        displayed_vector = length * mode[index] / maximum
        arrows.append({
            "id": f"display-{index + 1}-source-{source_index + 1}",
            "origin": (positions[index] - 0.5 * displayed_vector).tolist(),
            "origin_space": "cartesian",
            "vector": displayed_vector.tolist(),
            "direction_space": "cartesian",
            "color": ARROW_COLOR,
            "metadata": {
                "atom_index_1based": source_index + 1,
                "element": symbols[index],
                "relative_amplitude": float(relative[index]),
            },
        })
    group = {
        "id": "normal-mode",
        "name": "Normal-mode displacement",
        "magnitude_mode": "absolute",
        "viewport_policy": "include",
        "opacity": 1.0,
        "style": {
            "shaft_radius": 0.060,
            "head_length_ratio": 0.25,
            "head_radius_ratio": 2.3,
            "sides": 12,
            "flatshading": True,
        },
        "arrows": arrows,
    }
    return [group], {
        "selected_atom_indices_1based": [int(source_indices[index] + 1) for index in selected],
        "selected_atom_count": len(selected),
        "minimum_relative_amplitude": minimum_relative_amplitude,
        "maximum_display_arrow_length_ang": length,
    }


def font(size: int):
    path = Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf")
    return ImageFont.truetype(str(path), size) if path.is_file() else ImageFont.load_default()


def add_header_image(
    image: Image.Image,
    text: str,
    height: int = 54,
    font_size: int = 22,
) -> Image.Image:
    image = image.convert("RGB")
    canvas = Image.new("RGB", (image.width, image.height + height), "white")
    canvas.paste(image, (0, height))
    draw = ImageDraw.Draw(canvas)
    size = font_size
    selected_font = font(size)
    while draw.textbbox((0, 0), text, font=selected_font)[2] > image.width - 24 and size > 12:
        size -= 1
        selected_font = font(size)
    box = draw.textbbox((0, 0), text, font=selected_font)
    draw.text(((image.width - (box[2] - box[0])) / 2, 12), text, fill="#111111", font=selected_font)
    return canvas


def crop_frames(frames: list[Image.Image], padding: int) -> list[Image.Image]:
    x0, y0 = frames[0].width, frames[0].height
    x1 = y1 = 0
    for frame in frames:
        array = np.asarray(frame.convert("RGB"))
        ys, xs = np.where(np.any(array < 248, axis=2))
        if len(xs):
            x0, y0 = min(x0, int(xs.min())), min(y0, int(ys.min()))
            x1, y1 = max(x1, int(xs.max()) + 1), max(y1, int(ys.max()) + 1)
    if x1 <= x0 or y1 <= y0:
        raise ValueError("No visible content found in rendered frames")
    box = (
        max(0, x0 - padding), max(0, y0 - padding),
        min(frames[0].width, x1 + padding), min(frames[0].height, y1 + padding),
    )
    return [frame.crop(box) for frame in frames]


def pad_to_aspect(image: Image.Image, aspect: float) -> Image.Image:
    current = image.width / image.height
    if current > aspect:
        height = int(np.ceil(image.width / aspect))
        canvas = Image.new("RGB", (image.width, height), "white")
        canvas.paste(image, (0, (height - image.height) // 2))
    else:
        width = int(np.ceil(image.height * aspect))
        canvas = Image.new("RGB", (width, image.height), "white")
        canvas.paste(image, ((width - image.width) // 2, 0))
    return canvas


def render_frame(scene, style, vectors, camera, width, height, scale) -> Image.Image:
    figure = build_figure(scene, style, vector_overlays=vectors, include_interaction_traces=False)
    figure.layout.title = {"text": " "}
    figure.update_layout(
        scene_camera=camera,
        width=width,
        height=height,
        margin={"l": 0, "r": 0, "t": 6, "b": 6},
        paper_bgcolor="#FFFFFF",
    )
    image = figure.to_image(format="png", width=width, height=height, scale=scale)
    return Image.open(io.BytesIO(image)).convert("RGB")


def displaced_scene(base_scene: dict, positions: np.ndarray, displacement: np.ndarray) -> dict:
    scene = dict(base_scene)
    scene["draw_atoms"] = [dict(atom) for atom in base_scene["draw_atoms"]]
    for atom, position, delta in zip(scene["draw_atoms"], positions, displacement):
        atom["cart"] = np.asarray(position + delta)
    scene["bonds"] = [dict(bond) for bond in base_scene["bonds"]]
    for bond in scene["bonds"]:
        bond["start"] = np.asarray(scene["draw_atoms"][bond["i"]]["cart"])
        bond["end"] = np.asarray(scene["draw_atoms"][bond["j"]]["cart"])
    return scene


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path, help="*_gamma_modes.npz")
    parser.add_argument("--modes", type=int, nargs="+", required=True, help="One-based mode indices")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--label", help="Material label; defaults to archive stem")
    parser.add_argument("--assignments", nargs="*", help="Optional labels aligned with --modes")
    parser.add_argument("--camera-axis", choices=["a", "b", "c", "a*", "b*", "c*"], default="a")
    parser.add_argument(
        "--focus",
        choices=["all", "organic", "azide", "framework"],
        default="all",
        help="Render only the selected chemical component as an unwrapped cluster",
    )
    parser.add_argument("--arrow-length", type=float, default=2.0)
    parser.add_argument("--minimum-relative-amplitude", type=float, default=0.12)
    parser.add_argument("--heavy-atoms-only", action="store_true")
    parser.add_argument("--amplitude", type=float, default=0.35)
    parser.add_argument("--frames", type=int, default=16)
    parser.add_argument("--fps", type=float, default=10.0)
    parser.add_argument("--width", type=int, default=950)
    parser.add_argument("--height", type=int, default=760)
    parser.add_argument("--scale", type=int, default=2)
    parser.add_argument(
        "--png-scale",
        type=int,
        help="Static PNG render scale; defaults to --scale (use 6 for publication output)",
    )
    parser.add_argument("--png-dpi", type=int, default=600, help="DPI metadata embedded in static PNGs")
    parser.add_argument("--crop-padding", type=int, default=20)
    parser.add_argument("--aspect", type=float, default=1.60)
    parser.add_argument("--gif", action="store_true")
    parser.add_argument("--write-trajectories", action="store_true")
    args = parser.parse_args()
    if args.assignments and len(args.assignments) != len(args.modes):
        parser.error("--assignments must contain exactly one label per mode")
    if args.frames < 4 or args.amplitude <= 0 or args.fps <= 0:
        parser.error("--frames must be >= 4 and amplitude/fps must be positive")
    if args.arrow_length <= 0 or not 0 <= args.minimum_relative_amplitude <= 1:
        parser.error("arrow length must be positive and amplitude threshold must be in [0, 1]")
    if min(args.width, args.height, args.scale) <= 0 or args.crop_padding < 0 or args.aspect <= 0:
        parser.error("image dimensions, scale, and aspect must be positive; padding cannot be negative")
    png_scale = args.png_scale or args.scale
    if png_scale <= 0 or args.png_dpi <= 0:
        parser.error("PNG scale and DPI must be positive")

    data = load_archive(args.archive)
    atoms = build_atoms(data)
    group_names, participation = atom_groups(data)
    focus_indices = classified_source_indices(atoms)[args.focus]
    if args.focus == "all":
        payload, cell, matrix = payload_from_atoms(atoms)
        canonical_bond_pairs = None
        display_mode = "unit_cell"
        show_unit_cell = True
    else:
        payload, cell, matrix, canonical_bond_pairs = focused_payload(atoms, focus_indices)
        display_mode = "cluster"
        show_unit_cell = False
    camera = axis_camera(matrix, args.camera_axis)
    material = args.label or args.archive.name.removesuffix("_gamma_modes.npz")
    args.outdir.mkdir(parents=True, exist_ok=True)

    preset = {"style": {
        "show_title": False,
        "show_hydrogen": True,
        "show_labels": False,
        "show_axes": False,
        "show_unit_cell": show_unit_cell,
        "material": "mesh",
        "style": "ball_stick",
        "projection": "orthographic",
        "atom_scale": 0.72,
        "bond_radius": 0.105,
        "background": "#FFFFFF",
        "element_colors": ELEMENT_COLORS,
    }}
    scene = build_scene_from_atoms(
        name=material,
        title=material,
        atoms=payload,
        cell=cell,
        M=matrix,
        R=np.eye(3),
        show_hydrogen=True,
        preset=preset,
        display_mode=display_mode,
        ops=scene_ops(),
        bond_scale=1.15,
        bond_thresholds={
            ("Cu", "N"): 2.65 / 1.15,
            ("Pb", "N"): 3.10 / 1.15,
        },
        canonical_bond_pairs=canonical_bond_pairs,
        include_boundary_replicas=False,
    )
    positions = np.asarray([atom["cart"] for atom in scene["draw_atoms"]], dtype=float)
    symbols = [str(atom["elem"]) for atom in scene["draw_atoms"]]
    source_indices = np.asarray(
        [int(atom["_source_index"]) for atom in scene["draw_atoms"]], dtype=int
    )
    style = scene_style(scene, preset["style"])
    style["camera"] = camera
    summaries = []

    for item_index, mode_number in enumerate(args.modes):
        mode_index = mode_number - 1
        if not 0 <= mode_index < len(data["frequencies_thz"]):
            raise ValueError(f"Mode {mode_number} is outside the archive")
        source_mode = np.asarray(data["cartesian_displacements"][mode_index], dtype=float)
        mode = source_mode[source_indices]
        vectors, arrow_summary = mode_arrows(
            positions,
            symbols,
            mode,
            length=args.arrow_length,
            minimum_relative_amplitude=args.minimum_relative_amplitude,
            heavy_atoms_only=args.heavy_atoms_only,
            source_indices=source_indices,
        )
        frequency = float(data["frequencies_thz"][mode_index])
        dominant_group = max(group_names, key=lambda name: participation[name][mode_index])
        assignment = args.assignments[item_index] if args.assignments else dominant_group.replace("_", " ")
        header = f"{material} Γ mode {mode_number} · {frequency:.2f} THz · {assignment}"
        phases = np.linspace(
            0.0,
            2.0 * np.pi,
            args.frames if args.gif else 1,
            endpoint=False,
        )
        frames = [
            render_frame(
                displaced_scene(scene, positions, args.amplitude * np.sin(phase) * mode),
                style, vectors, camera, args.width, args.height, args.scale,
            )
            for phase in phases
        ]
        frames = [pad_to_aspect(add_header_image(frame, header), args.aspect) for frame in crop_frames(frames, args.crop_padding)]
        axis_slug = args.camera_axis.replace("*", "star")
        focus_slug = "" if args.focus == "all" else f"_{args.focus}"
        stem = (
            f"{material.replace(' ', '_')}_mode_{mode_number:03d}_"
            f"{frequency:.2f}THz_{axis_slug}{focus_slug}_mattervis"
        )
        png_path = args.outdir / f"{stem}.png"
        gif_path = args.outdir / f"{stem}.gif"
        if png_scale == args.scale:
            png_frame = frames[0]
        else:
            high_resolution = render_frame(
                displaced_scene(scene, positions, np.zeros_like(mode)),
                style, vectors, camera, args.width, args.height, png_scale,
            )
            ratio = png_scale / args.scale
            high_resolution = crop_frames(
                [high_resolution],
                max(1, round(args.crop_padding * ratio)),
            )[0]
            png_frame = pad_to_aspect(
                add_header_image(
                    high_resolution,
                    header,
                    height=max(1, round(54 * ratio)),
                    font_size=max(12, round(22 * ratio)),
                ),
                args.aspect,
            )
        png_frame.save(png_path, dpi=(args.png_dpi, args.png_dpi))
        if args.gif:
            frames[0].save(
                gif_path,
                save_all=True,
                append_images=frames[1:],
                duration=int(round(1000.0 / args.fps)),
                loop=0,
                disposal=2,
                optimize=False,
            )
        trajectory_path = None
        if args.write_trajectories:
            trajectory_path = args.outdir / f"{stem}.extxyz"
            trajectory = []
            trajectory_phases = np.linspace(0.0, 2.0 * np.pi, args.frames, endpoint=False)
            for phase in trajectory_phases:
                frame = atoms.copy()
                frame.positions = atoms.positions + args.amplitude * np.sin(phase) * source_mode
                frame.info.update({"mode_index_1based": mode_number, "frequency_thz": frequency})
                trajectory.append(frame)
            ase_write(trajectory_path, trajectory)
        summaries.append({
            "mode_index_1based": mode_number,
            "frequency_thz": frequency,
            "assignment_label": assignment,
            "dominant_group": dominant_group,
            "participation": {name: float(participation[name][mode_index]) for name in group_names},
            "png": str(png_path),
            "gif": str(gif_path) if args.gif else None,
            "trajectory": str(trajectory_path) if trajectory_path else None,
            "camera_axis": args.camera_axis,
            "focus": args.focus,
            "amplitude_angstrom": args.amplitude,
            "frame_count": args.frames,
            "fps": args.fps,
            "png_scale": png_scale,
            "png_dpi": args.png_dpi,
            "png_pixel_size": list(png_frame.size),
            **arrow_summary,
        })
        print(f"Saved {png_path}")
        if args.gif:
            print(f"Saved {gif_path}")

    axis_slug = args.camera_axis.replace("*", "star")
    mode_slug = "-".join(str(mode) for mode in args.modes)
    metadata_path = args.outdir / (
        f"{material.replace(' ', '_')}_gamma_modes_{mode_slug}_{axis_slug}.json"
    )
    metadata_path.write_text(json.dumps({
        "source_archive": str(args.archive),
        "mattervis_checkout": str(MATTERVIS),
        "material": material,
        "focus": args.focus,
        "note": "Assignment labels are user-provided or participation-based; visually inspect modes before mechanistic labeling.",
        "modes": summaries,
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Saved {metadata_path}")


if __name__ == "__main__":
    main()
