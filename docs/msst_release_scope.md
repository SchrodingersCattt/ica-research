# MSST release scope

`workflows/msst/in.lmp` is the sanitized LAMMPS/MSST input distributed with
this review bundle. It is an inspectable starting point, not a self-contained
job: it expects a user-supplied LAMMPS data file (`conf.lmp`), a compatible
DeepMD model, and a LAMMPS build with the MSST and DeepMD pair styles.

The template preserves the released review input, but its single-model and
time-integration settings should not be read as the final manuscript run
manifest. In particular, it uses a 72 A/ps shock velocity and one 0.0002 ps
timestep for both equilibration and shock stages; the manuscript describes
68 A/ps, a 1 fs (0.001 ps) NPT stage followed by a 0.1 fs (0.0001 ps) MSST
stage, and a 20 ps equilibration. The template equilibrates for 1,000 steps
and writes coordinates every 100 steps, whereas the manuscript states 1,000
steps for coordinate output. These values need to be reconciled against the
archived production run before claiming exact MSST reproduction.

The public AIS Square artifact documented in `models/README.md` is the frozen
DPA3 model used for phonon and elastic inference. The manuscript's MSST methods
describe a committee of compressed `se_atten_v2` models; those private or
cluster-specific model files are not included in this Git repository.

Reaction-network analysis is described in `docs/pseudocode.md` and requires an
independently installed ReacNetGenerator with the bond-perception and atom
origin rules stated in the manuscript. The release does not include raw MSST
trajectories or a ReacNetGenerator executable.

Accordingly, this repository supports code inspection and reproduction of the
portable examples and analysis workflows. Exact reproduction of the reported
MSST trajectory and reaction-network figures additionally requires the
external inputs, model committee, software versions, and run manifest.
