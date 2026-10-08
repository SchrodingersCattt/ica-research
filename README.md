# ICA Research computational workflows

This review bundle contains the source code and small input/output examples used for the computational parts of the CA-series manuscript. It covers DPA3/DeepMD fine-tuning and inference, finite-displacement phonons, elastic constants, LAMMPS MSST input, and CP2K/NCI post-processing.

This is a local review release (`0.1.0-review`). It has not been pushed to GitHub and does not contain the large external programs themselves. Multiwfn and VMD are cited as external dependencies only.

## Quick map

| Area | Entry points | External engine |
| --- | --- | --- |
| Phonons | `workflows/phonon/calc_phonon.py`, `calc_phonon_tetrahedron.py` | DeepMD-kit, ASE, Phonopy, SeeK-path |
| Elasticity | `workflows/elastic/ase_elastic.py` | DeepMD-kit, ASE, Pymatgen |
| Shock simulation | `workflows/msst/in.lmp` | LAMMPS MSST, DeepMD pair style |
| Electron density/NCI | `workflows/dft_nci/cp2k_neT4_CA.inp`, plotting scripts | CP2K, Multiwfn, VMD/Tachyon |
| Fine-tuning | `workflows/finetune/` | DeepMD-kit and a site-specific launcher supplied locally |

## System requirements

The recorded production environment was Linux with Python 3.11, CUDA-capable GPUs, and a DeepMD environment. The exact observed Python package versions are in `docs/software_versions.md`. A CPU-only machine can run the plotting and stored-result analysis examples, but the full phonon, elastic, MSST, and fine-tuning workflows require the corresponding external engines and model files.

Do not put cluster paths, credentials, `.aissq/`, or site-specific launcher settings in this directory. Copy `config.local.example.yml` to a local ignored file and fill in paths on the machine where a calculation will run.

## Installation and run-time estimates

For the Python analysis-only examples, a normal desktop installation is expected to take 5–15 minutes. The stored elastic summary smoke test should finish in under one minute. The phonon example is a GPU calculation and can take minutes to hours depending on the supercell and model. Full fine-tuning and MSST runs are cluster jobs and are not run by the smoke test.

Create an environment from `environment.yml` where compatible wheels are available, then run:

```bash
python -m compileall workflows tools
python tools/smoke_test.py
```

The smoke test validates imports that are installed, reads the example structures and elastic JSON files, checks the sanitized training template, and performs static checks on CP2K/LAMMPS inputs. It deliberately does not submit jobs or launch external GUI software.

## Reproduction workflows

The complete command-level description is in `docs/reproduce.md`; algorithm summaries and the manuscript mapping are in `docs/pseudocode.md`. The phonon mode notes are retained in `workflows/phonon/README_modes.md`.

Typical phonon invocation, after providing a compatible model and environment:

```bash
python workflows/phonon/calc_phonon.py \
  --systems CA LA NEt4-CA \
  --model "$ICA_MODEL_PATH" \
  --cif-dir examples/structures \
  --output-dir output/phonon
```

Elastic analysis can use stored JSON results without rerunning the expensive calculation:

```bash
python workflows/elastic/ase_elastic.py --mode plot --output-dir output/elastic
```

Fine-tuning configuration is documented but its launcher is intentionally site-neutral. The original cloud/DLC identifiers and absolute NAS paths were removed from this bundle.

## Models and external software

`models/README.md` records the frozen DPA3 inference model provenance and checksum. The model itself is kept outside Git until the review is complete; later publication can use the aissq artifact workflow. Full checkpoints and training datasets are likewise external.

`docs/software_versions.md` lists the versions supported by evidence in the source logs. CP2K is referenced rather than compiled; the manuscript's CP2K version discrepancy is called out there for final confirmation. Multiwfn and VMD/Tachyon are referenced only and are not redistributed.

The exact public-client capability record is in `docs/aissq-explorer.md`; the authenticated AIS Square upload sequence and artifact record are in `docs/aissquare-upload-api.md`.

## License

Original scripts in this bundle are released under the MIT License. Third-party programs, libraries, model weights, and data remain under their own licenses; see `THIRD_PARTY_NOTICES.md`.
