# MSST and reaction-network workflow

`in.lmp` is the LAMMPS/MSST input. `conf.lmp` is a public, reproducible
starting configuration generated from `examples/structures/NEt4-CA.cif` by a
4 × 3 × 3 replication. It is not the private post-NPT production snapshot.
The file contains no velocities; `in.lmp` initializes them before the run.

After installing ReacNetGenerator, run the sanitized wrapper on a locally
stored trajectory:

```bash
python workflows/msst/run_reacnet.py \
  --trajectory /path/to/dump_msst_freq-1000.lammpstrj \
  --output-dir output/reacnet
```

The wrapper uses ASE-based bond perception, excludes Cu–N connectivity with
`Cu-N:0`, disables hidden-Markov filtering, and records a path-free manifest.
Use `--dry-run` to inspect the command without reading a trajectory. The
trajectory and ReacNetGenerator executable remain external dependencies.

The manuscript's source-resolved N₂ fractions additionally require assigning
each N₂ molecule from the timestep-zero nitrogen identities (azide,
organic-cation, or mixed). That provenance step is documented in the methods
and should be applied to the local ReacNet sidecars; no trajectory is bundled
with this repository.
