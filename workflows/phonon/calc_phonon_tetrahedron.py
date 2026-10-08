import itertools
import math
import os
import json
from ase import io
from ase.io import read
from deepmd.calculator import DP as DeepMD
from ase.optimize import (
    LBFGS, FIRE, BFGS, BFGSLineSearch
)
from ase.filters import ExpCellFilter, FrechetCellFilter 
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
import seekpath
import numpy as np
from ase import Atoms
import glob
from tqdm import tqdm
import logging
import random

np.random.seed(42)
random.seed(42)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', filename='phonon_calculation_tetrahedron.log')

# Function to read the structure from a CIF file
def read_structure(file_path: str):
    """
    Read a structure from a CIF file using ASE.
    """
    return io.read(file_path)

# Function to relax the structure using DeepMD calculator
def relax_structure(atoms, model_path: str):
    """
    Fully relax the structure using the DeepMD neural network potential.
    """
    atoms.calc = DeepMD(model=model_path)
    ecf = FrechetCellFilter(atoms)
    opt = BFGS(ecf)
    opt.run(fmax=0.005)  # Convergence threshold for forces
    return atoms

# Function to find high-symmetry points using seekpath
def get_high_symmetry_points(atoms):
    """
    Use seekpath to find the high-symmetry points and path in the Brillouin zone.
    """
    cell = atoms.get_cell()
    positions = atoms.get_scaled_positions()
    numbers = atoms.get_atomic_numbers()
    structure = (cell, positions, numbers)

    kpath_data = seekpath.get_path(structure)
    return kpath_data, structure

def calculate_phonon_spectrum(atoms, supercell_matrix, model_path, displacement_distance):
    # Convert ASE Atoms to PhonopyAtoms
    ph_atoms = PhonopyAtoms(
        symbols=atoms.get_chemical_symbols(),
        cell=atoms.get_cell(),
        scaled_positions=atoms.get_scaled_positions()
    )
    # Create Phonopy instance using the PhonopyAtoms object and the supercell matrix
    phonon = Phonopy(
        ph_atoms, 
        supercell_matrix=supercell_matrix,
        # symprec=1e-3
    )
    phonon.generate_displacements(distance=displacement_distance, is_diagonal=False)
    
    _calc = DeepMD(model=model_path)
    logging.info("Calculating forces")
    forcesets = []

    for frame in tqdm(phonon.supercells_with_displacements):
        frame_atom = Atoms(
            cell=frame.cell,
            symbols=frame.symbols,
            scaled_positions=frame.scaled_positions,
            pbc=True,
        )
        frame_atom.calc = _calc
        forces = frame_atom.get_forces()
        drift_force = forces.sum(axis=0)
        for force in forces:
            force -= drift_force / forces.shape[0]
        forcesets.append(forces)

    # Use the forces attribute instead of the deprecated set_forces method
    phonon.forces = forcesets
    phonon.produce_force_constants()
    phonon.symmetrize_force_constants()

    return phonon

from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections

def save_phonon_band_structure_plot(phonon, kpath_data, structure, save_path, nx, ny, nz):
    path = [list(kpath_data['point_coords'].values())]
    labels = list(kpath_data['point_coords'].keys())
    
    logging.info("Generating band structure")
    phonon.auto_band_structure(
        npoints=101,
        write_yaml=True,
        filename=f"{save_path.split('.png')[0]}_band.yaml"
    )
    plt_band = phonon.plot_band_structure()
    plt_band.savefig(f"{save_path.split('.png')[0]}_band.png", dpi=300)
    logging.info("Generating band structure with DOS")
    
    # MAX_ATOMS_RATIO = 0.5 
    MAX_ATOMS_RATIO = 1.0
    MIN_ATOM = 1 

    pdos_indices = [
        random.sample(
            [i for i, val in enumerate(structure[2]) if val == value],
            k=max(MIN_ATOM, math.ceil(MAX_ATOMS_RATIO * len([i for i, val in enumerate(structure[2]) if val == value])))
        )
        for value in set(structure[2])
    ]
    logging.info(f"{save_path.split('/')[-1].split('.')[0]}: {structure[2]}\n{pdos_indices}")
    logging.info(f"nx, ny, nz = {nx}, {ny}, {nz}")
    
    # Using a denser mesh for more accurate tetrahedron DOS calculation
    phonon.run_mesh(
        [
            nx * 4, 
            ny * 4, 
            nz * 4
        ], 
        with_eigenvectors=True, 
        is_mesh_symmetry=False,
    )
    
    # Calculate DOS using tetrahedron method
    phonon.run_total_dos(use_tetrahedron_method=True, sigma=0.01)
    phonon.plot_band_structure_and_dos().savefig(save_path.replace('.png', '_with_tetrahedron_tdos.png'), dpi=300)
    
    # Calculate projected DOS using tetrahedron method
    phonon.run_projected_dos(use_tetrahedron_method=True, sigma=0.01)
    phonon.plot_band_structure_and_dos(
        pdos_indices=pdos_indices
    ).savefig(save_path.replace('.png', '_with_tetrahedron_pdos.png'), dpi=300)
    
    # Write DOS data to files with tetrahedron-specific names
    phonon.write_total_dos(
        filename=f"{save_path.split('.png')[0]}_tetrahedron_tdos.dat"
    )
    phonon.write_projected_dos(
        filename=f"{save_path.split('.png')[0]}_tetrahedron_pdos.dat"
    )
    
    phonon.run_thermal_properties(temperatures=300)
    thermal = phonon.get_thermal_properties_dict()
    with open(save_path.replace('.png', '_tetrahedron_thermal.json'), 'w') as f:
        json.dump(thermal, f, default=lambda o: o.tolist() if isinstance(o, np.ndarray) else o)
    os.system(
        f"phonopy-bandplot --gnuplot {save_path.split('.png')[0]}_band.yaml > {save_path.split('.png')[0]}_tetrahedron_band.dat"
    )


# Main function to orchestrate the workflow
def calc_phonon(system: str, nx: int, ny: int, nz: int, cif_dir: str, model_path: str, output_dir: str):
    """
    Main workflow to calculate the phonon spectrum for a given system.
    """
    
    supercell_matrix = [[nx, 0, 0], [0, ny, 0], [0, 0, nz]] 
    os.makedirs(output_dir, exist_ok=True)
    cif_path = os.path.join(cif_dir, f"{system}.cif")
    _atoms = read_structure(cif_path)

    # Relax the structure
    logging.info(f"Relaxing structure for {system}...")
    atoms = relax_structure(_atoms, model_path)

    # Find high-symmetry points
    logging.info("Finding high-symmetry points...")
    kpath_data, structure = get_high_symmetry_points(atoms)
    logging.info(f"High-symmetry path: {kpath_data['path']}")

    # Set up phonon calculation
    logging.info("Calculating phonon spectrum...")
    displacement_distance = 0.01
    phonon = calculate_phonon_spectrum(atoms, supercell_matrix, model_path, displacement_distance)

    # Save the phonon band structure plot
    plot_path = os.path.join(output_dir, f"{system}_tetrahedron.png")
    save_phonon_band_structure_plot(phonon, kpath_data, structure, plot_path, nx, ny, nz)

    logging.info(f"Phonon spectrum saved as {plot_path}")

# Entry point
if __name__ == "__main__":
    confs = glob.glob("confs/**.cif")
    for conf in confs:
        atom = read(conf)
        a, b, c, _, _, _ = atom.cell.cellpar()
        system = conf.split("/")[-1].split(".")[0]
        # nx, ny, nz = 3, 12, 5
        N = 50
        nx, ny, nz = (
            max(1, N // a), 
            max(1, N // b), 
            max(1, N // c)
        )
        logging.info(f"{system} {a} {b} {c}; {nx} {ny} {nz}")
        calc_phonon(
            system, 
            nx, ny, nz, 
            model_path="dpa3.pth", 
            cif_dir='confs', 
            output_dir="output"
        )