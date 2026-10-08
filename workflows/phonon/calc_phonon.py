import itertools
import math
import os
import json
import argparse
import hashlib
from datetime import datetime, timezone
from ase import io
from ase.io import read
from ase.neighborlist import neighbor_list, natural_cutoffs
from deepmd.calculator import DP as DeepMD
from ase.optimize import (
    LBFGS, FIRE, BFGS, BFGSLineSearch
)
from ase.filters import ExpCellFilter, FrechetCellFilter 
from phonopy import Phonopy
from phonopy import __version__ as phonopy_version
from phonopy.file_IO import write_force_constants_to_hdf5
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

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', filename='phonon_calculation.log')


def file_sha256(path: str):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def classify_atom_groups(atoms):
    """Classify atoms into framework, organic-cation, and azide groups.

    C and H atoms are organic. Nitrogen bonded to carbon under periodic
    minimum-image geometry is organic; all remaining nitrogen is azide.
    This avoids atom-order-specific lists while retaining a transparent
    chemistry-based definition for NEt4-CA.
    """
    symbols = np.asarray(atoms.get_chemical_symbols())
    cutoffs = natural_cutoffs(atoms, mult=1.15)
    left, right = neighbor_list("ij", atoms, cutoffs)
    organic_n = set()
    nitrogen_adjacency = {
        index: set() for index, symbol in enumerate(symbols) if symbol == "N"
    }
    for i, j in zip(left, right):
        if {symbols[i], symbols[j]} == {"C", "N"}:
            organic_n.add(i if symbols[i] == "N" else j)
        if symbols[i] == symbols[j] == "N":
            nitrogen_adjacency[i].add(j)

    organic = [
        i for i, symbol in enumerate(symbols)
        if symbol in {"C", "H"} or i in organic_n
    ]
    azide = [i for i, symbol in enumerate(symbols) if symbol == "N" and i not in organic_n]
    remaining = set(azide)
    components = []
    while remaining:
        root = remaining.pop()
        component = {root}
        stack = [root]
        while stack:
            current = stack.pop()
            for neighbour in nitrogen_adjacency[current] & remaining:
                remaining.remove(neighbour)
                component.add(neighbour)
                stack.append(neighbour)
        components.append(component)
    if any(len(component) != 3 for component in components):
        sizes = sorted(len(component) for component in components)
        raise ValueError(f"Non-organic nitrogen does not form N3 components: {sizes}")
    assigned = set(organic) | set(azide)
    framework = [i for i in range(len(atoms)) if i not in assigned]
    groups = {"framework": framework, "azide": azide}
    if organic:
        groups["organic_cation"] = organic
    if sorted(index for values in groups.values() for index in values) != list(range(len(atoms))):
        raise ValueError("Atom-group classification is not a complete, disjoint partition")
    return groups


def phase_aligned_modes(eigenvector_matrix, masses):
    """Return mode-major eigenvectors and real Cartesian display modes."""
    mode_count = eigenvector_matrix.shape[1]
    mass_weighted = np.empty((mode_count, len(masses), 3), dtype=np.complex128)
    cartesian = np.empty((mode_count, len(masses), 3), dtype=float)
    imaginary_fraction = np.empty(mode_count, dtype=float)

    for mode_index in range(mode_count):
        vector = np.asarray(eigenvector_matrix[:, mode_index], dtype=np.complex128)
        pivot = int(np.argmax(np.abs(vector)))
        if abs(vector[pivot]) > 0:
            vector *= np.exp(-1j * np.angle(vector[pivot]))
        shaped = vector.reshape(-1, 3)
        mass_weighted[mode_index] = shaped
        imaginary_fraction[mode_index] = np.linalg.norm(shaped.imag) / max(np.linalg.norm(shaped), 1e-30)
        display = shaped.real / np.sqrt(masses[:, None])
        maximum = np.linalg.norm(display, axis=1).max()
        cartesian[mode_index] = display / maximum if maximum > 0 else display
    return mass_weighted, cartesian, imaginary_fraction


def save_gamma_mode_data(
    phonon,
    atoms,
    system,
    output_dir,
    model_path,
    supercell_matrix,
    displacement_distance,
):
    """Persist Gamma eigenvectors, force constants, groups, and provenance."""
    phonon.run_qpoints([[0.0, 0.0, 0.0]], with_eigenvectors=True)
    qpoints = phonon.get_qpoints_dict()
    frequencies = np.asarray(qpoints["frequencies"][0], dtype=float)
    eigenvector_matrix = np.asarray(qpoints["eigenvectors"][0], dtype=np.complex128)
    masses = np.asarray(phonon.primitive.masses, dtype=float)
    primitive_symbols = list(phonon.primitive.symbols)
    primitive_positions = np.asarray(phonon.primitive.scaled_positions, dtype=float)
    ase_positions = np.asarray(atoms.get_scaled_positions(), dtype=float)
    wrapped_delta = primitive_positions - ase_positions
    wrapped_delta -= np.rint(wrapped_delta)
    if primitive_symbols != atoms.get_chemical_symbols() or not np.allclose(wrapped_delta, 0.0, atol=1e-7):
        raise ValueError("Phonopy primitive atom order does not match the relaxed ASE structure")
    if eigenvector_matrix.shape != (3 * len(atoms), 3 * len(atoms)):
        raise ValueError(f"Unexpected Gamma eigenvector shape: {eigenvector_matrix.shape}")
    mass_weighted, cartesian, imaginary_fraction = phase_aligned_modes(eigenvector_matrix, masses)
    groups = classify_atom_groups(atoms)
    mass_weighted_participation = {
        name: np.sum(np.abs(mass_weighted[:, indices, :]) ** 2, axis=(1, 2))
        for name, indices in groups.items()
    }
    mass_weighted_total = np.sum(list(mass_weighted_participation.values()), axis=0)
    mass_weighted_participation = {
        name: values / mass_weighted_total for name, values in mass_weighted_participation.items()
    }
    cartesian_participation = {
        name: np.sum(cartesian[:, indices, :] ** 2, axis=(1, 2))
        for name, indices in groups.items()
    }
    cartesian_total = np.sum(list(cartesian_participation.values()), axis=0)
    cartesian_participation = {
        name: values / cartesian_total for name, values in cartesian_participation.items()
    }

    prefix = os.path.join(output_dir, system)
    io.write(f"{prefix}_relaxed.cif", atoms)
    write_force_constants_to_hdf5(
        phonon.force_constants,
        filename=f"{prefix}_force_constants.hdf5",
        physical_unit="eV/angstrom^2",
    )
    np.savez_compressed(
        f"{prefix}_gamma_modes.npz",
        frequencies_thz=frequencies,
        mass_weighted_eigenvectors=mass_weighted,
        cartesian_displacements=cartesian,
        imaginary_fraction=imaginary_fraction,
        symbols=np.asarray(atoms.get_chemical_symbols()),
        masses=masses,
        cell=np.asarray(atoms.cell),
        scaled_positions=np.asarray(atoms.get_scaled_positions()),
        group_names=np.asarray(list(groups)),
        participation=np.column_stack([mass_weighted_participation[name] for name in groups]),
        mass_weighted_participation=np.column_stack([
            mass_weighted_participation[name] for name in groups
        ]),
        cartesian_participation=np.column_stack([
            cartesian_participation[name] for name in groups
        ]),
    )
    metadata = {
        "system": system,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "mode_indexing": "one-based in reports; zero-based in NumPy arrays",
        "frequency_unit": "THz",
        "displacement_unit": "dimensionless; each mode normalized to max atom norm = 1",
        "eigenvector_convention": "Phonopy mass-weighted eigenvectors, phase aligned by largest component",
        "cartesian_conversion": "real(phase-aligned eigenvector) / sqrt(atomic mass)",
        "participation_metrics": {
            "mass_weighted_participation": "sum(abs(Phonopy eigenvector)^2) by atom group",
            "cartesian_participation": "sum(abs(mass-unweighted Cartesian displacement)^2) by atom group",
            "participation": "backward-compatible alias of mass_weighted_participation",
        },
        "model_path": os.path.abspath(model_path),
        "model_sha256": file_sha256(model_path),
        "phonopy_version": phonopy_version,
        "supercell_matrix": np.asarray(supercell_matrix).tolist(),
        "displacement_distance_angstrom": displacement_distance,
        "atom_groups_1based": {
            name: [index + 1 for index in indices] for name, indices in groups.items()
        },
        "max_gamma_imaginary_fraction": float(imaginary_fraction.max()),
    }
    with open(f"{prefix}_gamma_modes.json", "w", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2)
        stream.write("\n")
    logging.info("Saved Gamma modes and force constants for %s", system)

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

def save_band_structure_only(phonon, save_path):
    """Persist and plot the band path without running the expensive DOS mesh."""
    logging.info("Generating band structure only")
    phonon.auto_band_structure(
        npoints=101,
        write_yaml=True,
        filename=f"{save_path.split('.png')[0]}_band.yaml"
    )
    plt_band = phonon.plot_band_structure()
    plt_band.savefig(f"{save_path.split('.png')[0]}_band.png", dpi=300)
    os.system(
        f"phonopy-bandplot --gnuplot {save_path.split('.png')[0]}_band.yaml > "
        f"{save_path.split('.png')[0]}_band.dat"
    )


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
    phonon.run_mesh(
        [
            nx * 10, 
            ny * 10,  
            nz * 10
        ], 
        # shift=[0.25, 0.25, 0.25],
        with_eigenvectors=True, 
        is_mesh_symmetry=False,
        # is_gamma_center=True,
    )
    phonon.run_total_dos(sigma=0.01)
    phonon.plot_band_structure_and_dos().savefig(save_path.replace('.png', '_with_tdos.png'), dpi=300)
    
    phonon.run_projected_dos(sigma=0.01)
    phonon.plot_band_structure_and_dos(
        pdos_indices=pdos_indices
    ).savefig(save_path.replace('.png', '_with_pdos.png'), dpi=300)
    phonon.write_total_dos(
        filename=f"{save_path.split('.png')[0]}_tdos.dat"
    )
    phonon.write_projected_dos(
        filename=f"{save_path.split('.png')[0]}_pdos.dat"
    )
    phonon.run_thermal_properties(temperatures=300)
    thermal = phonon.get_thermal_properties_dict()
    with open(save_path.replace('.png', '_thermal.json'), 'w') as f:
        json.dump(thermal, f, default=lambda o: o.tolist() if isinstance(o, np.ndarray) else o)
    os.system(
        f"phonopy-bandplot --gnuplot {save_path.split('.png')[0]}_band.yaml > {save_path.split('.png')[0]}_band.dat"
    )


# Main function to orchestrate the workflow
def calc_phonon(
    system: str,
    nx: int,
    ny: int,
    nz: int,
    cif_dir: str,
    model_path: str,
    output_dir: str,
    gamma_only: bool = False,
    skip_relaxation: bool = False,
    band_only: bool = False,
):
    """
    Main workflow to calculate the phonon spectrum for a given system.
    """
    
    supercell_matrix = [[nx, 0, 0], [0, ny, 0], [0, 0, nz]] 
    os.makedirs(output_dir, exist_ok=True)
    cif_path = os.path.join(cif_dir, f"{system}.cif")
    _atoms = read_structure(cif_path)

    # Relax the structure
    if skip_relaxation:
        logging.info("Using input structure without relaxation for %s", system)
        atoms = _atoms
    else:
        logging.info(f"Relaxing structure for {system}...")
        atoms = relax_structure(_atoms, model_path)

    # Find high-symmetry points
    logging.info("Finding high-symmetry points...")
    kpath_data, structure = get_high_symmetry_points(atoms)
    logging.info(f"High-symmetry path: {kpath_data['path']}")

    # Set up phonon calculation
    logging.info("Calculating phonon spectrum...")
    # if system == 'Et':
    #     displacement_distance = 0.03
    # else:
    #     displacement_distance = 0.01
    displacement_distance = 0.01
    phonon = calculate_phonon_spectrum(atoms, supercell_matrix, model_path, displacement_distance)

    save_gamma_mode_data(
        phonon,
        atoms,
        system,
        output_dir,
        model_path,
        supercell_matrix,
        displacement_distance,
    )
    if gamma_only:
        logging.info("Gamma-only mode requested; skipping band and DOS calculations for %s", system)
        return

    # Save the phonon band structure plot
    plot_path = os.path.join(output_dir, f"{system}.png")
    if band_only:
        save_band_structure_only(phonon, plot_path)
    else:
        save_phonon_band_structure_plot(phonon, kpath_data, structure, plot_path, nx, ny, nz)

    logging.info(f"Phonon spectrum saved as {plot_path}")

# Entry point
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate phonons and persist Gamma eigenmodes")
    parser.add_argument(
        "--systems",
        nargs="+",
        choices=["CA", "LA", "NEt4-CA"],
        help="Only process these systems (default: every CIF in confs/)",
    )
    parser.add_argument("--model", default="dpa3.pth")
    parser.add_argument("--cif-dir", default="confs")
    parser.add_argument("--output-dir", default="output")
    parser.add_argument(
        "--gamma-only",
        action="store_true",
        help="Persist Gamma modes and force constants, then skip band and DOS calculations",
    )
    parser.add_argument(
        "--skip-relaxation",
        action="store_true",
        help="Use CIF geometries directly; intended for retry-relaxed structures extracted from band YAML",
    )
    parser.add_argument(
        "--band-only",
        action="store_true",
        help="Generate Gamma modes and the full band path, but skip DOS and thermal calculations",
    )
    args = parser.parse_args()
    confs = glob.glob(os.path.join(args.cif_dir, "*.cif"))
    for conf in confs:
        atom = read(conf)
        a, b, c, _, _, _ = atom.cell.cellpar()
        system = conf.split("/")[-1].split(".")[0]
        if args.systems and system not in args.systems:
            continue
        # nx, ny, nz = 3, 12, 5
        N = 40
        nx, ny, nz = (
            max(1, N // a), 
            max(1, N // b), 
            max(1, N // c)
        )
        logging.info(f"{system} {a} {b} {c}; {nx} {ny} {nz}")
        calc_phonon(
            system, 
            nx, ny, nz, 
            model_path=args.model,
            cif_dir=args.cif_dir,
            output_dir=args.output_dir,
            gamma_only=args.gamma_only,
            skip_relaxation=args.skip_relaxation,
            band_only=args.band_only,
        )
