# Computational pseudocode and manuscript mapping

## Phonons and thermal properties

1. Read a CIF with ASE and attach the DPA3/DeepMD calculator.
2. Relax the cell and atomic positions unless the input is an already relaxed structure.
3. Construct a Phonopy supercell and finite displacements.
4. Evaluate forces for every displaced supercell with the ML potential.
5. Build force constants, obtain the SeeK-path high-symmetry path, and calculate bands and projected DOS.
6. Export native band/PDOS/thermal files and plot-ready tables.

This corresponds to the manuscript paragraph beginning “Phonon calculations were performed using the finite-displacement supercell method”.

## Elastic constants

1. Relax the reference structure with the same ML potential.
2. Generate six independent normal/shear strain states at four amplitudes.
3. Relax atoms at fixed cell shape and collect the six stress components.
4. Fit stress versus strain, assemble the 6x6 tensor, and save the raw fit points.
5. Compute Voigt, Reuss, Hill, Young, Poisson and universal-anisotropy quantities.

This corresponds to the manuscript paragraph beginning “The full second-order elastic constant tensor”.

## Fine-tuning and validation

1. Load a DeepMD JSON template.
2. Substitute data roots, model dimensions, loss prefactors, learning-rate schedule and batch size from local configuration.
3. Run a site-specific launcher only after a dry-run has printed the resolved command.
4. Read logs and validation systems, evaluate energy/force/virial predictions, and report RMSE values.

## MSST and reaction analysis

1. Build and equilibrate the replicated crystal with LAMMPS and the DeepMD pair style.
2. Apply the MSST shock along the specified crystallographic direction.
3. Save thermodynamic quantities and atomic trajectories.
4. Run ReacNetGenerator externally, with the documented bond-perception and atom-origin rules, and export reaction summaries.

## DFT/NCI

1. Run the supplied CP2K input externally to obtain periodic electron-density cube files.
2. Use Multiwfn externally to calculate IRI and sign(lambda2)rho fields.
3. Use VMD/Tachyon externally to render the isosurfaces. The executables are not part of this repository.
