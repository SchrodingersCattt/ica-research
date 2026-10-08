# Nature Research software checklist audit

This table maps the supplied Code and Software Submission Checklist to the
current repository. “Partial” means the repository states the limitation
explicitly and identifies the evidence still needed before claiming full
reproducibility.

| Checklist item | Status | Evidence in this repository |
| --- | --- | --- |
| Compiled software and/or source code | Pass | `workflows/` contains the source scripts and calculation inputs; large external programs are referenced, not redistributed. |
| Small demo dataset | Pass | `examples/structures/`, `examples/*_elastic.json`, and `examples/phonon_expected/`. |
| README system requirements | Pass | `README.md` records Linux, Python 3.11, NVIDIA A800 80 GB, driver 535.129.03 and CUDA 12.2. |
| All dependencies and versions | Partial | Exact ASE, DeepMD-kit and Phonopy versions are recorded; Pymatgen, SeeK-path, dpdata, LAMMPS, ReacNetGenerator and VMD/Tachyon still need production evidence. |
| Non-standard hardware | Pass | GPU model, driver and CUDA versions are recorded in `README.md` and `docs/software_versions.md`. |
| Installation guide | Pass | `environment.yml` and the installation commands in `README.md` are provided; DeepMD-kit is called out as hardware-specific. |
| Typical installation time | Pass | `README.md` gives the 5–15 minute analysis-only estimate. |
| Demo instructions and run time | Pass | `tools/smoke_test.py` and the expected runtime are documented in `README.md`. |
| Instructions for use | Pass | `docs/reproduce.md` and workflow READMEs provide commands and inputs. |
| Reproduction instructions | Pass | `docs/reproduce.md` and `docs/pseudocode.md`. |
| Software license | Pass | `LICENSE` and `THIRD_PARTY_NOTICES.md`. |
| Open source repository link | Pass | `README.md` links to `SchrodingersCattt/ica-research`. |
| Detailed functionality / pseudocode | Pass | `docs/pseudocode.md` maps each workflow to the manuscript methods. |

## Remaining release gates

1. Record exact production versions for Pymatgen, SeeK-path, dpdata, LAMMPS,
   ReacNetGenerator, VMD and Tachyon.
2. Resolve the manuscript CP2K 2023.1 wording against the job module recorded
   as CP2K 2023.2.
3. Ask an unfamiliar colleague to install the environment and run the smoke
   test, as recommended by the supplied checklist.
