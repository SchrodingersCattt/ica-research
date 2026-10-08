# Reproduction guide

1. Create an environment from `environment.yml` and install a compatible DeepMD-kit build.
2. Obtain the frozen DPA3 model from the artifact URL in `models/README.md`.
3. Copy `config.local.example.yml` to ignored `config.local.yml` and set local paths.
4. Run `python tools/smoke_test.py` before any expensive job.
5. Run the phonon and elastic workflows with the example structures; compare generated native files with the expected-output directory.
6. Submit fine-tuning, CP2K and LAMMPS jobs only through the local site configuration. The public scripts contain no cluster identifiers.
7. Perform Multiwfn, VMD/Tachyon and ReacNetGenerator steps with independently installed software, recording the exact versions in the local run manifest.
