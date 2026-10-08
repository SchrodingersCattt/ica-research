# Software and version evidence

The production phonon log records the following exact Python runtime versions:

```text
Python 3.11
ASE 3.28.0
DeepMD-kit 3.1.3.dev21+gb98f6c596
Phonopy 2.41.1
NVIDIA driver 535.129.03
CUDA 12.2
GPU NVIDIA A800-SXM4 80 GB
```

These are the versions used for the stored phonon results. The conda file
contains the portable Python dependencies; DeepMD-kit remains a separate
hardware-specific installation because its CUDA, TensorFlow and compiler
compatibility must match the execution host.

| Software | Role | Version recorded | Evidence/status |
| --- | --- | --- | --- |
| Python | runtime | 3.11 | recorded from the source production log; the raw cluster log is not redistributed |
| ASE | atomic structures and calculators | 3.28.0 | recorded from the source production log; the raw cluster log is not redistributed |
| DeepMD-kit | DPA3 model inference/training | 3.1.3.dev21+gb98f6c596 | source production log; hardware/build specific |
| Phonopy | phonons and thermal properties | 2.41.1 | recorded from the source production log; the raw cluster log is not redistributed |
| Pymatgen | elasticity analysis | version not recorded | record from the production environment before release |
| SeeK-path | high-symmetry paths | version not recorded | record from the production environment before release |
| dpdata | DeepMD data conversion | version not recorded | record from the production environment before release |
| CP2K | periodic electron density | 2023.2-gcc-11.4.0-openmpi-5.0.0-ch4 module | exact module string is in `workflows/dft_nci/cp2k.slurm.example`; reconcile against manuscript's 2023.1 wording |
| LAMMPS | MSST molecular dynamics | not recorded | cite upstream release; do not invent a version |
| ReacNetGenerator | reaction network analysis | not recorded | cite upstream release; do not invent a version |
| Multiwfn | IRI/NCI scalar fields | 3.8 in manuscript methods | external software reference only |
| VMD/Tachyon | NCI/IRI rendering | not recorded | external software reference only |

No external executable is compiled or redistributed by this repository. The final public release must add immutable upstream/release links and fill entries currently marked “not recorded” when evidence is available.
