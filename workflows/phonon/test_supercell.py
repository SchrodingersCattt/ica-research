import os
import math
import json
import numpy as np
import logging
import random
from tqdm import tqdm
from ase.build import bulk
from ase import Atoms
from ase.optimize import BFGS
from ase.filters import FrechetCellFilter
from deepmd.calculator import DP as DeepMD
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
import seekpath

np.random.seed(42)
random.seed(42)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
def relax_structure(atoms, model_path):
    atoms.calc = DeepMD(model=model_path)
    ecf = FrechetCellFilter(atoms)
    opt = BFGS(ecf)
    opt.run(fmax=0.01)
    return atoms

def get_high_symmetry_points(atoms):
    cell = atoms.get_cell()
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()
    structure = (cell, positions, numbers)
    kpath_data = seekpath.get_path(structure)
    return kpath_data, structure

def calculate_phonon_spectrum(atoms, supercell_matrix, model_path, displacement_distance):
    ph_atoms = PhonopyAtoms(
        symbols=atoms.get_chemical_symbols(),
        cell=atoms.get_cell(),
        scaled_positions=atoms.get_scaled_positions()
    )
    phonon = Phonopy(ph_atoms, supercell_matrix=supercell_matrix)
    phonon.generate_displacements(distance=displacement_distance, is_diagonal=False)

    calc = DeepMD(model=model_path)
    logging.info("Calculating forces")
    forcesets = []

    for frame in tqdm(phonon.supercells_with_displacements):
        displaced = Atoms(
            cell=frame.cell,
            symbols=frame.symbols,
            scaled_positions=frame.scaled_positions,
            pbc=True,
        )
        displaced.calc = calc
        forces = displaced.get_forces()
        drift = forces.sum(axis=0)
        forces -= drift / forces.shape[0]
        forcesets.append(forces)

    phonon.forces = forcesets
    phonon.produce_force_constants()
    phonon.symmetrize_force_constants()

    return phonon

def save_phonon_outputs(phonon, kpath_data, structure, save_dir, nx, ny, nz):
    import matplotlib.pyplot as plt

    path = [list(kpath_data['point_coords'].values())]
    phonon.auto_band_structure(npoints=101)

    # DOS
    phonon.run_mesh([nx * 4, ny * 4, nz * 4], with_eigenvectors=True, is_mesh_symmetry=False)
    phonon.run_total_dos()

    # Projected DOS
    pdos_indices = [
        random.sample(
            [i for i, val in enumerate(structure[2]) if val == z],
            k=max(1, math.ceil(0.5 * len([i for i, val in enumerate(structure[2]) if val == z])))
        )
        for z in set(structure[2])
    ]
    phonon.run_projected_dos()
    # Thermal properties
    phonon.run_thermal_properties(temperatures=range(100, 1001, 100))
    thermal = phonon.get_thermal_properties_dict()
    logging.info(thermal)
def run_phonon_for_fcc_cu(nx, ny, nz, model_path):
    logging.info(f"Running for FCC Cu with supercell {nx}x{ny}x{nz}")
    from ase.io import read
    atoms = read("confs/CA.cif")  # typical lattice constant
    atoms.set_pbc(True)
    atoms = relax_structure(atoms, model_path)
    kpath_data, structure = get_high_symmetry_points(atoms)
    phonon = calculate_phonon_spectrum(atoms, [[nx, 0, 0], [0, ny, 0], [0, 0, nz]], model_path, displacement_distance=0.01)
    save_dir = f"output/Cu_{nx}x{ny}x{nz}"
    save_phonon_outputs(phonon, kpath_data, structure, save_dir, nx, ny, nz)

if __name__ == "__main__":
    model_path = "dpa3.pth"
    supercells = [(2, 2, 2), (3, 3, 3)]
    for nx, ny, nz in supercells:
        run_phonon_for_fcc_cu(nx, ny, nz, model_path)
