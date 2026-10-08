# Γ-point vibration-mode workflow

> `phonon_retry_2609` is an isolated copy. The original `phonon_retry` tree is
> not used as a write target. Its immutable comparison files are copied under
> `reference_retry/`.

## 2609 reproduction and deliverables

`confs_retry_relaxed/` contains structures extracted directly from the original
retry `*_band.yaml` files. This avoids an additional relaxation and allows the
new force-constant calculation to reproduce the retry spectral features as
closely as possible.

Run:

```bash
python calc_phonon.py --systems CA LA NEt4-CA \
  --model "$ICA_MODEL_PATH" --cif-dir ../../examples/structures \
  --output-dir output_2609 --skip-relaxation --band-only
```

The calculation writes only to:

- `output_2609/`: regenerated full band paths, Γ eigenvectors, force constants,
  relaxed-geometry copies, provenance, and candidate-mode tables;
- `comparison/`: side-by-side retry/2609 spectra and numerical metrics;
- `mode_media/`: selected PNG/GIF/EXTXYZ vibration deliverables.
- `deliverable_600dpi/`: publication-resolution PNGs (600-DPI metadata and
  roughly 4,000–6,000 pixels across), plus matching lightweight GIFs,
  trajectories, and comparison metrics.

After the calculation, inspect `output_2609/gamma_mode_candidates.csv`, select
the clearest organic and azide modes, and render them with
`render_gamma_modes.py`. Do not assign torsion/libration labels until the GIFs
have been visually inspected.

## PDOS grouping and canonical figures

The canonical combined outputs are written directly as:

- `band_and_dos.png`
- `Azides_phonon_and_dos.pdf`
- `vib_pdos.png`

No `corrected` or legacy-backup suffix is used. Regenerating a figure replaces
the same canonical filename.

For ICA, Cu, organic-cation, and azide PDOS groups are noncontiguous in the
underlying atom order. `plot_band_and_vdos.py` and `plot_vdos.py` must therefore
sum the explicit zero-based atom-index lists after adding one for the frequency
column in `*_pdos.dat`. Grouping contiguous blocks merely from each list length
is incorrect: it mixes azide N into the organic curve and C/H into the azide
curve. That error produces a spurious organic contribution near 61–64 THz.

With the explicit-index grouping, the ICA assignments are:

- 35–46 THz: predominantly organic-cation modes (integrated contribution
  82.12% organic and 17.88% azide in the current retry PDOS);
- 60–66 THz: azide-stretching manifold (98.15% azide and 1.85% organic);
- 61–64 THz: likewise azide dominated (98.12% azide).

Thus the high-frequency feature near 61–64 THz must be labelled as N₃⁻
stretching, not as an Et₄N⁺ organic mode.

This workflow generates defensible mode-resolved PNG/GIF media for CA, LA, and
NEt₄–CA (ICA), using the same `dpa3.pth` model and supercell policy as
`calc_phonon.py`.

## 1. Recalculate and persist eigenmodes

Run from this directory with the project environment:

```bash
python calc_phonon.py \
  --systems CA LA NEt4-CA --gamma-only
```

This is the expensive step: it repeats relaxation and finite-displacement force
calculations. `--gamma-only` skips the subsequent band-path and dense-mesh DOS
stages. It writes per material:

- `output/<system>_relaxed.cif`
- `output/<system>_force_constants.hdf5`
- `output/<system>_gamma_modes.npz`
- `output/<system>_gamma_modes.json`

The NPZ contains Γ frequencies, phase-aligned mass-weighted eigenvectors,
mass-unweighted Cartesian display displacements, both mass-weighted and
Cartesian-motion atom-group participation, and structure data. The JSON records
model hash and calculation provenance.

## 2. Rank candidate modes

```bash
python analyze_gamma_modes.py \
  output/CA_gamma_modes.npz \
  output/LA_gamma_modes.npz \
  output/NEt4-CA_gamma_modes.npz
```

Inspect `output/gamma_mode_candidates.csv` and
`output/gamma_mode_candidates.json`. Participation identifies candidate modes;
it does not by itself justify labels such as torsion, libration, deformation,
or stretching.

Recommended initial inspection windows:

- ICA low-frequency organic: 0.1–12 THz
- ICA intermediate organic: 35–46 THz, especially near 41.67 THz
- Azide stretch controls: 60–70 THz in all three materials

## 3. Render selected modes

Example only—replace mode indices after inspecting the ranking table:

```bash
python render_gamma_modes.py \
  output/NEt4-CA_gamma_modes.npz \
  --modes 541 619 \
  --label ICA-1 \
  --assignments "organic deformation" "azide stretch" \
  --outdir mode_media/ICA-1 \
  --camera-axis a --gif --write-trajectories
```

The renderer uses an optional local MatterVis checkout and writes PNG, optional GIF,
optional EXTXYZ trajectories, and JSON metadata. Render promising modes from a
second camera axis before choosing the publication still. Camera axes are
included in filenames, so alternate views do not overwrite each other. Set
`CA_SERIES_MATTERVIS_ROOT` if the reference checkout is moved.

For high-resolution static figures, add `--png-scale 6 --png-dpi 600`. The
ordinary `--scale` continues to control GIF frames, preventing unnecessarily
large animations while keeping the static PNG suitable for publication.

## Publication checks

1. Confirm regenerated Γ frequencies agree with the corresponding first Γ
   entries in `output/*_band.yaml`.
2. Confirm group participation sums to one for every mode.
3. Visually inspect all assignment labels. Degenerate modes can rotate within a
   subspace and must not be over-interpreted individually.
4. Use one consistent data set (`phonon_retry`) for the VDOS and mode media.
5. Do not report the existing 39.97, 38.23, and 62.66 THz onset triplet until a
   single documented threshold reproduces all three values.
