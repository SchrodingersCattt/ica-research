# Software and version evidence

| Software | Role | Version recorded | Evidence/status |
| --- | --- | --- | --- |
| ASE | atomic structures and calculators | 3.28.0 | `data/phonon_retry_2609/run_2609.log` |
| DeepMD-kit | DPA3 model inference/training | 3.1.3.dev21+gb98f6c596 | `data/phonon_retry_2609/run_2609.log`; hardware/build specific |
| Phonopy | phonons and thermal properties | 2.41.1 | `data/phonon_retry_2609/run_2609.log` |
| Pymatgen | elasticity analysis | version not recorded | record from the production environment before release |
| SeeK-path | high-symmetry paths | version not recorded | record from the production environment before release |
| dpdata | DeepMD data conversion | version not recorded | record from the production environment before release |
| CP2K | periodic electron density | 2023.2 in SLURM module template | reconcile against the manuscript's 2023.1 wording using job evidence |
| LAMMPS | MSST molecular dynamics | not recorded | cite upstream release; do not invent a version |
| ReacNetGenerator | reaction network analysis | not recorded | cite upstream release; do not invent a version |
| Multiwfn | IRI/NCI scalar fields | 3.8 in manuscript methods | external software reference only |
| VMD/Tachyon | NCI/IRI rendering | not recorded | external software reference only |

No external executable is compiled or redistributed by this repository. The final public release must add immutable upstream/release links and fill entries currently marked “not recorded” when evidence is available.
