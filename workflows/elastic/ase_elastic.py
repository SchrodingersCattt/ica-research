"""
Elastic constant calculation and visualization.

Usage
-----
# Step 1: run ML calculations and save results to JSON
python ase_elastic.py --mode calc [--cif-glob ../phonon/confs/**cif] [--output-dir output]

# Step 2: plot from saved JSON (tweak freely without re-running MD)
python ase_elastic.py --mode plot [--output-dir output]

# Run both in one go
python ase_elastic.py --mode all
"""

import argparse
import csv
import json
import logging
import os
from pathlib import Path
from typing import List, TypedDict, Optional

import glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

EV_A3_TO_GPA = 160.21766208  # eV/Å³ to GPa

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    filename="elastic.log",
)


# ---------------------------------------------------------------------------
# Type
# ---------------------------------------------------------------------------

class ElasticResult(TypedDict):
    material: str
    elastic_tensor: List[List[float]]   # 6×6, GPa (Voigt)
    bulk_modulus: float                  # Voigt B, GPa
    shear_modulus: float                 # Voigt G, GPa
    youngs_modulus: float                # GPa
    poisson_ratio: float
    # raw fit data for stress-strain plots: shape [6][n_points]
    fit_strains: List[List[float]]       # fit_strains[i] = strains for strain state i
    fit_stresses: List[List[List[float]]] # fit_stresses[i][j] = stress component j for strain state i
    fit_r2: List[List[float]]            # fit_r2[i][j] = R² for Cij fit


# ---------------------------------------------------------------------------
# Calculation helpers  (only imported when --mode calc/all)
# ---------------------------------------------------------------------------

def _get_elastic_tensor_from_strains(strains, stresses, eq_stress=None, tol=1e-7):
    """
    Returns (elastic_tensor, fit_data) where fit_data contains raw strains,
    stresses, slopes and R² for each (i,j) pair.
    """
    from pymatgen.analysis.elasticity import ElasticTensor
    from pymatgen.analysis.elasticity.elastic import get_strain_state_dict

    strain_states = [tuple(ss) for ss in np.eye(6)]
    ss_dict = get_strain_state_dict(
        strains, stresses,
        eq_stress=eq_stress,
        add_eq=(eq_stress is not None),
    )
    c_ij = np.zeros((6, 6))
    fit_strains  = []   # [6] list of 1D arrays
    fit_stresses = []   # [6] list of (n, 6) arrays
    fit_r2       = []   # [6][6]

    for ii in range(6):
        strain_arr = ss_dict[strain_states[ii]]["strains"]  # (n, 6)
        stress_arr = ss_dict[strain_states[ii]]["stresses"]  # (n, 6)
        fit_strains.append(strain_arr[:, ii].tolist())
        fit_stresses.append(stress_arr.tolist())
        r2_row = []
        for jj in range(6):
            fit = np.polyfit(strain_arr[:, ii], stress_arr[:, jj], 1, full=True)
            slope = fit[0][0]
            c_ij[ii, jj] = slope
            y_true = stress_arr[:, jj]
            y_pred = slope * strain_arr[:, ii]
            ss_res = np.sum((y_true - y_pred) ** 2)
            ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
            r2 = 1 - ss_res / ss_tot if ss_tot != 0 else float("nan")
            r2_row.append(float(r2))
            print(f"  Fit σ[{jj}] vs ε[{ii}]: slope={slope:.4f}, R²={r2:.4f}")
        fit_r2.append(r2_row)

    elastic_tensor = ElasticTensor.from_voigt(c_ij)
    print("Elastic tensor (GPa):")
    print(elastic_tensor.voigt * EV_A3_TO_GPA)
    return elastic_tensor.zeroed(tol), fit_strains, fit_stresses, fit_r2


def calculate_elastic_constants(
    cif_file: Path,
    model_path: Path,
    head: str = "Omat24",
    norm_strains=None,
    shear_strains=None,
    fmax: float = 0.01,
) -> ElasticResult:
    from ase.constraints import ExpCellFilter
    from ase.io import read, write
    from ase.optimize import FIRE
    from deepmd.calculator import DP
    from pymatgen.analysis.elasticity import DeformedStructureSet, Strain, Stress
    from pymatgen.io.ase import AseAtomsAdaptor

    if norm_strains is None:
        norm_strains = np.linspace(-0.01, 0.01, 4)
    if shear_strains is None:
        shear_strains = np.linspace(-0.01, 0.01, 4)

    material = Path(cif_file).stem
    print(f"\n{'='*60}\nProcessing {cif_file}\n{'='*60}")

    initial_atoms = read(str(cif_file))
    calc = DP(model=str(model_path), head=head)
    initial_atoms.calc = calc

    ecf = ExpCellFilter(initial_atoms, mask=[True]*6)
    FIRE(ecf, logfile="relax.log").run(fmax=fmax)
    relaxed_atoms = ecf.atoms
    write("relax.cif", relaxed_atoms)

    structure = AseAtomsAdaptor.get_structure(relaxed_atoms)
    deformed_set = DeformedStructureSet(
        structure, norm_strains=norm_strains, shear_strains=shear_strains
    )

    stresses = []
    for ds in deformed_set:
        atoms = ds.to_ase_atoms()
        atoms.calc = calc
        FIRE(atoms, logfile=None).run(fmax=fmax)
        stresses.append(atoms.get_stress(voigt=False) * EV_A3_TO_GPA)

    strains = [Strain.from_deformation(d) for d in deformed_set.deformations]

    relaxed_atoms.calc = calc
    eq_stress = Stress(relaxed_atoms.get_stress(voigt=False) * EV_A3_TO_GPA)

    et, fit_strains, fit_stresses, fit_r2 = _get_elastic_tensor_from_strains(
        strains, stresses, eq_stress=eq_stress
    )

    B = float(et.k_voigt)
    G = float(et.g_voigt)
    E = float(9 * B * G / (3 * B + G))
    nu = float(0.5 * (3 * B - 2 * G) / (3 * B + G))

    return {
        "material": material,
        "elastic_tensor": [[float(et.voigt[i][j]) for j in range(6)] for i in range(6)],
        "bulk_modulus": B,
        "shear_modulus": G,
        "youngs_modulus": E,
        "poisson_ratio": nu,
        "fit_strains":   fit_strains,
        "fit_stresses":  fit_stresses,
        "fit_r2":        fit_r2,
    }


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def save_result(result: ElasticResult, output_dir: str) -> None:
    path = os.path.join(output_dir, f"{result['material']}_elastic.json")
    with open(path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"Saved → {path}")


def load_results(output_dir: str) -> List[ElasticResult]:
    paths = sorted(glob.glob(os.path.join(output_dir, "*_elastic.json")))
    if not paths:
        raise FileNotFoundError(
            f"No *_elastic.json found in {output_dir}. Run --mode calc first."
        )
    results = []
    for p in paths:
        with open(p) as f:
            results.append(json.load(f))
    print(f"Loaded {len(results)} result(s) from {output_dir}")
    return results


# ---------------------------------------------------------------------------
# Elasticity analysis: VRH, anisotropy, soft-direction statistics
# ---------------------------------------------------------------------------

def _voigt_reuss_hill(C: np.ndarray) -> dict:
    S = np.linalg.inv(C)
    B_V = (C[0,0]+C[1,1]+C[2,2] + 2*(C[0,1]+C[0,2]+C[1,2])) / 9.0
    G_V = (C[0,0]+C[1,1]+C[2,2] - C[0,1]-C[0,2]-C[1,2]
           + 3*(C[3,3]+C[4,4]+C[5,5])) / 15.0
    B_R = 1.0 / (S[0,0]+S[1,1]+S[2,2] + 2*(S[0,1]+S[0,2]+S[1,2]))
    G_R = 15.0 / (4*(S[0,0]+S[1,1]+S[2,2])
                  - 4*(S[0,1]+S[0,2]+S[1,2])
                  + 3*(S[3,3]+S[4,4]+S[5,5]))
    B_H = (B_V + B_R) / 2.0
    G_H = (G_V + G_R) / 2.0
    A_U = 5 * G_V / G_R + B_V / B_R - 6.0
    return dict(B_V=B_V, B_R=B_R, B_H=B_H,
                G_V=G_V, G_R=G_R, G_H=G_H,
                delta_G=G_V-G_R, A_U=A_U)


def _compliance_4th(S: np.ndarray):
    voigt = {
        (0,0):0,(1,1):1,(2,2):2,
        (1,2):3,(2,1):3,(0,2):4,(2,0):4,(0,1):5,(1,0):5,
    }
    def s(i,j,k,l):
        p,q = voigt[(i,j)], voigt[(k,l)]
        f = 1.0
        if p > 2: f *= 0.5
        if q > 2: f *= 0.5
        return S[p,q] * f
    return s


def _gmin_gmax_exact(S: np.ndarray, n: np.ndarray):
    """Exact G_min and G_max via 2×2 eigenvalue problem."""
    if abs(n[2]) < 0.9:
        e1 = np.cross(n, [0,0,1])
    else:
        e1 = np.cross(n, [1,0,0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1); e2 /= np.linalg.norm(e2)
    basis = [e1, e2]
    s = _compliance_4th(S)
    A = np.zeros((2,2))
    for a in range(2):
        for b in range(2):
            val = 0.0
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        for l in range(3):
                            val += s(i,j,k,l)*n[i]*basis[a][j]*n[k]*basis[b][l]
            A[a,b] = 4.0 * val
    eigvals = np.linalg.eigvalsh(A)
    G_min = 1.0/eigvals[-1] if eigvals[-1] > 1e-12 else float("nan")
    G_max = 1.0/eigvals[0]  if eigvals[0]  > 1e-12 else float("nan")
    return G_min, G_max


def _build_shear_surface_exact(S: np.ndarray, n_theta: int = 80, n_phi: int = 160):
    """Returns (X,Y,Z, G_min_arr, G_max_arr, n_vectors) with r=1/G_min."""
    theta = np.linspace(0, np.pi, n_theta)
    phi   = np.linspace(0, 2*np.pi, n_phi)
    THETA, PHI = np.meshgrid(theta, phi)
    Nx = np.sin(THETA)*np.cos(PHI)
    Ny = np.sin(THETA)*np.sin(PHI)
    Nz = np.cos(THETA)
    G_min_arr = np.zeros_like(THETA)
    G_max_arr = np.zeros_like(THETA)
    for i in range(n_phi):
        for j in range(n_theta):
            n_vec = np.array([Nx[i,j], Ny[i,j], Nz[i,j]])
            gmin, gmax = _gmin_gmax_exact(S, n_vec)
            G_min_arr[i,j] = gmin
            G_max_arr[i,j] = gmax
    R = 1.0 / G_min_arr
    return R*Nx, R*Ny, R*Nz, G_min_arr, G_max_arr, Nx, Ny, Nz


def compute_soft_stats(G_min_arr: np.ndarray, Nx, Ny, Nz,
                       G_c: Optional[float] = None):
    """
    Compute:
      - G_min_global: global minimum G_min
      - n_star: direction of global minimum (unit vector)
      - f_soft: fraction of directions with G_min < G_c
      - G_std: std of G_min distribution
    G_c defaults to the 20th percentile of all G_min values across all materials
    (set externally for cross-material comparison).
    """
    flat = G_min_arr.flatten()
    flat = flat[np.isfinite(flat)]
    G_min_global = float(np.min(flat))
    idx = np.argmin(G_min_arr)
    i_min, j_min = np.unravel_index(idx, G_min_arr.shape)
    n_star = np.array([Nx[i_min, j_min], Ny[i_min, j_min], Nz[i_min, j_min]])

    if G_c is None:
        G_c = float(np.percentile(flat, 20))

    f_soft = float(np.sum(flat < G_c) / len(flat))
    G_std  = float(np.std(flat))

    return dict(
        G_min_global=G_min_global,
        n_star=n_star.tolist(),
        G_c=G_c,
        f_soft=f_soft,
        G_std=G_std,
    )


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

COLORS = ["#4C72B0", "#DD8452", "#55A868"]


def plot_stress_strain_fits(result: ElasticResult, output_dir: str) -> None:
    """
    6-panel figure: one panel per strain state (ε_ii or ε_ij).
    Each panel shows the diagonal stress component σ_jj vs ε_ii scatter + fit line,
    with R² annotated. Only the diagonal Cii fit is shown (most informative).
    """
    material = result["material"]
    if "fit_strains" not in result or result["fit_strains"] is None:
        print(f"  No fit data for {material}, skipping stress-strain plot.")
        return

    C = np.array(result["elastic_tensor"])
    voigt_labels = ["11", "22", "33", "23", "13", "12"]
    strain_labels = [f"ε$_{{{'11' if i==0 else '22' if i==1 else '33' if i==2 else '23' if i==3 else '13' if i==4 else '12'}}}$"
                     for i in range(6)]

    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    axes = axes.flatten()

    for ii in range(6):
        ax = axes[ii]
        eps = np.array(result["fit_strains"][ii])
        sig_all = np.array(result["fit_stresses"][ii])  # (n, 6)
        r2_row  = result["fit_r2"][ii]

        off_diag_labeled = False  # only add legend entry once for off-diagonal
        for jj in range(6):
            sig_j = sig_all[:, jj]
            slope = C[ii, jj]
            eps_line = np.linspace(eps.min(), eps.max(), 100)
            sig_line = slope * eps_line
            r2 = r2_row[jj]
            if jj == ii:
                # Diagonal Cii: highlighted in red
                ax.scatter(eps, sig_j, color="#C44E52", s=40, zorder=5,
                           label=f"σ$_{{{voigt_labels[jj]}}}$ data")
                ax.plot(eps_line, sig_line, color="#C44E52", lw=2.0, zorder=4,
                        label=f"fit: $C_{{{voigt_labels[ii]}{voigt_labels[jj]}}}$={slope:.1f} GPa\n$R^2$={r2:.4f}")
            else:
                # Off-diagonal Cij coupling terms: gray, one legend entry
                lbl = "coupling $C_{ij}$ ($j\\neq i$)" if not off_diag_labeled else None
                ax.scatter(eps, sig_j, color="gray", s=15, alpha=0.4, zorder=2,
                           label=lbl)
                ax.plot(eps_line, sig_line, color="gray", lw=0.8, alpha=0.4, zorder=1)
                off_diag_labeled = True

        ax.axhline(0, color="black", lw=0.5, ls="--")
        ax.axvline(0, color="black", lw=0.5, ls="--")
        ax.set_xlabel(f"Strain {strain_labels[ii]}", fontsize=10)
        ax.set_ylabel("Stress (GPa)", fontsize=10)
        ax.set_title(f"Strain state {ii+1}: {strain_labels[ii]}", fontsize=10)
        ax.legend(fontsize=7.5, loc="upper left")
        ax.grid(True, alpha=0.3)

    fig.suptitle(f"{material} — Stress–Strain Linear Fits", fontsize=13, y=1.01)
    fig.tight_layout()
    out = os.path.join(output_dir, f"{material}_stress_strain_fits.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved stress-strain fits → {out}")


def plot_shear_3d(result: ElasticResult, output_dir: str) -> dict:
    """
    3D shear compliance surface r=1/G_min(n).
    Returns soft_stats dict.
    """
    from mpl_toolkits.mplot3d import Axes3D  # noqa

    material = result["material"]
    C = np.array(result["elastic_tensor"])

    eigvals_C = np.linalg.eigvalsh(C)
    if np.any(eigvals_C <= 0):
        print(f"  WARNING: {material} not positive definite (min eig={eigvals_C.min():.4f}). Skip 3D.")
        return {}

    S = np.linalg.inv(C)
    print(f"  Computing exact 3D shear surface for {material}...")
    X, Y, Z, G_min_arr, G_max_arr, Nx, Ny, Nz = _build_shear_surface_exact(
        S, n_theta=80, n_phi=160
    )

    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection="3d")
    G_norm = (G_min_arr - G_min_arr.min()) / (G_min_arr.max() - G_min_arr.min() + 1e-12)
    ax.plot_surface(X, Y, Z, facecolors=plt.cm.coolwarm(G_norm),
                    rstride=1, cstride=1, alpha=0.92, linewidth=0, antialiased=True)

    sm = plt.cm.ScalarMappable(cmap="coolwarm",
                                norm=plt.Normalize(G_min_arr.min(), G_min_arr.max()))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.55, pad=0.1)
    cbar.set_label("$G_{\\min}$ (GPa)", fontsize=10)

    ax.set_xlabel("X (GPa$^{-1}$·scale)", fontsize=8)
    ax.set_ylabel("Y (GPa$^{-1}$·scale)", fontsize=8)
    ax.set_zlabel("Z (GPa$^{-1}$·scale)", fontsize=8)
    ax.set_title(
        f"{material}\nShear compliance surface $1/G_{{\\min}}(\\mathbf{{n}})$\n"
        f"$G_{{\\min}}$: {G_min_arr.min():.2f}–{G_min_arr.max():.2f} GPa",
        fontsize=10
    )
    ax.set_box_aspect([1,1,1])
    fig.tight_layout()
    out = os.path.join(output_dir, f"{material}_shear_3d.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved 3D shear surface → {out}")

    return G_min_arr, G_max_arr, Nx, Ny, Nz


def plot_soft_direction_summary(results: List[ElasticResult],
                                all_soft_stats: list,
                                output_dir: str) -> None:
    """
    Two-panel figure:
    Left:  G_min_global bar chart (absolute softest shear)
    Right: f_soft bar chart (fraction of soft directions, same G_c for all)
    Both panels annotate the softest direction n* as approximate Miller index.
    """
    materials = [r["material"] for r in results]
    G_min_g = [s["G_min_global"] for s in all_soft_stats]
    f_soft  = [s["f_soft"]       for s in all_soft_stats]
    G_c     = all_soft_stats[0]["G_c"]  # same threshold for all
    n_stars = [s["n_star"]       for s in all_soft_stats]
    G_stds  = [s["G_std"]        for s in all_soft_stats]

    def _approx_miller(n):
        """Approximate Miller index from unit vector."""
        n = np.array(n)
        n = n / np.abs(n).max()
        # round to nearest 0.5 then scale
        rounded = np.round(n * 2) / 2
        # find smallest integer representation
        from math import gcd
        from functools import reduce
        nums = [int(round(x * 2)) for x in rounded]
        g = reduce(gcd, [abs(x) for x in nums if x != 0])
        if g == 0: g = 1
        hkl = [x // g for x in nums]
        return f"[{hkl[0]} {hkl[1]} {hkl[2]}]"

    n = len(materials)
    x = np.arange(n)
    w = 0.5

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(max(10, 3.5*n), 5))

    # ── Left: G_min_global ─────────────────────────────────────────────────
    bars1 = ax1.bar(x, G_min_g, w, color=COLORS[:n], edgecolor="white", alpha=0.9)
    for bar, v, ns, gs in zip(bars1, G_min_g, n_stars, G_stds):
        miller = _approx_miller(ns)
        ax1.text(bar.get_x() + bar.get_width()/2, v + 0.05,
                 f"{v:.2f} GPa\n{miller}\nσ={gs:.2f}",
                 ha="center", va="bottom", fontsize=8)

    ax1.set_xticks(x); ax1.set_xticklabels(materials, fontsize=11)
    ax1.set_ylabel("$G_{\\min}^{\\rm global}$ (GPa)", fontsize=11)
    ax1.set_title("Global Minimum Shear Modulus", fontsize=12)
    ax1.grid(axis="y", linestyle="--", alpha=0.4)
    ax1.set_axisbelow(True)

    # ── Right: f_soft ──────────────────────────────────────────────────────
    bars2 = ax2.bar(x, [v*100 for v in f_soft], w,
                    color=COLORS[:n], edgecolor="white", alpha=0.9)
    for bar, v in zip(bars2, f_soft):
        ax2.text(bar.get_x() + bar.get_width()/2, v*100 + 0.3,
                 f"{v*100:.1f}%", ha="center", va="bottom", fontsize=9)

    ax2.set_xticks(x); ax2.set_xticklabels(materials, fontsize=11)
    ax2.set_ylabel(f"$f_{{\\rm soft}}$ (%) — directions with $G_{{\\min}} < {G_c:.1f}$ GPa",
                   fontsize=10)
    ax2.set_title("Soft Shear Direction Fraction $f_{\\rm soft}$", fontsize=12)
    ax2.grid(axis="y", linestyle="--", alpha=0.4)
    ax2.set_axisbelow(True)

    fig.tight_layout()
    out = os.path.join(output_dir, "soft_direction_summary.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved soft direction summary → {out}")


def save_summary_csv(results: List[ElasticResult], vrh_list: list,
                     all_soft_stats: list, output_dir: str) -> None:
    path = os.path.join(output_dir, "elastic_summary.csv")
    fields = ["material",
              "B_V", "G_V", "G_R", "G_H", "delta_G",
              "E", "nu", "pugh_BG", "A_U", "A_G",
              "G_min_global", "G_std", "f_soft", "G_c", "n_star"]
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for r, v, s in zip(results, vrh_list, all_soft_stats):
            writer.writerow({
                "material":    r["material"],
                "B_V":         f"{r['bulk_modulus']:.3f}",
                "G_V":         f"{v['G_V']:.3f}",
                "G_R":         f"{v['G_R']:.3f}",
                "G_H":         f"{v['G_H']:.3f}",
                "delta_G":     f"{v['delta_G']:.3f}",
                "E":           f"{r['youngs_modulus']:.3f}",
                "nu":          f"{r['poisson_ratio']:.4f}",
                "pugh_BG":     f"{r['bulk_modulus']/v['G_H']:.4f}",
                "A_U":         f"{v['A_U']:.4f}",
                "A_G":         f"{v['G_V']/v['G_R']:.4f}",
                "G_min_global":f"{s['G_min_global']:.3f}",
                "G_std":       f"{s['G_std']:.3f}",
                "f_soft":      f"{s['f_soft']:.4f}",
                "G_c":         f"{s['G_c']:.3f}",
                "n_star":      str(s["n_star"]),
            })
    print(f"Saved summary CSV → {path}")


def print_summary_table(results: List[ElasticResult], vrh_list: list,
                        all_soft_stats: list) -> None:
    hdr = (f"{'Material':<12} {'G_V':>6} {'G_R':>6} {'G_H':>6} "
           f"{'G_min':>7} {'G_std':>6} {'f_soft':>7} {'B/G_H':>7} {'A^U':>7}")
    print(f"\n{'='*len(hdr)}")
    print(hdr)
    print(f"{'-'*len(hdr)}")
    for r, v, s in zip(results, vrh_list, all_soft_stats):
        print(
            f"{r['material']:<12} "
            f"{v['G_V']:>6.2f} "
            f"{v['G_R']:>6.2f} "
            f"{v['G_H']:>6.2f} "
            f"{s['G_min_global']:>7.3f} "
            f"{s['G_std']:>6.3f} "
            f"{s['f_soft']*100:>6.1f}% "
            f"{r['bulk_modulus']/v['G_H']:>7.4f} "
            f"{v['A_U']:>7.4f}"
        )
    print(f"{'='*len(hdr)}\n")


# ---------------------------------------------------------------------------
# Mode runners
# ---------------------------------------------------------------------------

def run_calc(args) -> None:
    model_path = Path(args.model)
    cif_files = sorted(glob.glob(args.cif_glob))
    if not cif_files:
        raise FileNotFoundError(f"No CIF files matched: {args.cif_glob}")
    os.makedirs(args.output_dir, exist_ok=True)

    for cif_file in cif_files:
        try:
            result = calculate_elastic_constants(
                Path(cif_file), model_path,
                head=args.head,
                fmax=args.fmax,
            )
            save_result(result, args.output_dir)
        except Exception as e:
            logging.error(f"Failed {cif_file}: {e}", exc_info=True)
            print(f"FAILED {cif_file}: {e}")


def run_plot(args) -> None:
    os.makedirs(args.output_dir, exist_ok=True)
    results = load_results(args.output_dir)

    # VRH for all materials
    vrh_list = [_voigt_reuss_hill(np.array(r["elastic_tensor"])) for r in results]

    # ── Per-material: stress-strain fits + 3D surface ──────────────────────
    all_G_min_arrs = []
    all_Nx, all_Ny, all_Nz = [], [], []

    for r in results:
        plot_stress_strain_fits(r, args.output_dir)
        ret = plot_shear_3d(r, args.output_dir)
        if ret:
            G_min_arr, G_max_arr, Nx, Ny, Nz = ret
            all_G_min_arrs.append(G_min_arr)
            all_Nx.append(Nx); all_Ny.append(Ny); all_Nz.append(Nz)
        else:
            all_G_min_arrs.append(None)
            all_Nx.append(None); all_Ny.append(None); all_Nz.append(None)

    # ── Determine shared G_c (20th percentile of all G_min values) ─────────
    all_flat = []
    for arr in all_G_min_arrs:
        if arr is not None:
            all_flat.extend(arr.flatten().tolist())
    all_flat = [v for v in all_flat if np.isfinite(v)]
    G_c_shared = float(np.percentile(all_flat, 20))
    print(f"\nShared G_c (20th percentile across all materials) = {G_c_shared:.3f} GPa")

    # ── Compute soft stats with shared G_c ─────────────────────────────────
    all_soft_stats = []
    for r, arr, Nx, Ny, Nz in zip(results, all_G_min_arrs, all_Nx, all_Ny, all_Nz):
        if arr is not None:
            stats = compute_soft_stats(arr, Nx, Ny, Nz, G_c=G_c_shared)
        else:
            stats = {"G_min_global": float("nan"), "n_star": [0,0,0],
                     "G_c": G_c_shared, "f_soft": float("nan"), "G_std": float("nan")}
        all_soft_stats.append(stats)

    # ── Multi-material summary plots ────────────────────────────────────────
    plot_soft_direction_summary(results, all_soft_stats, args.output_dir)

    # ── Data export ────────────────────────────────────────────────────────
    save_summary_csv(results, vrh_list, all_soft_stats, args.output_dir)
    print_summary_table(results, vrh_list, all_soft_stats)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Elastic constant calculation and visualization.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--mode", choices=["calc","plot","all"], default="all")
    p.add_argument("--output-dir", default="output")
    p.add_argument("--cif-glob", default="../phonon/confs/**cif")
    p.add_argument(
        "--model",
        default=os.environ.get("ICA_MODEL_PATH", "dpa3.pth"),
    )
    p.add_argument("--head", default="Omat24")
    p.add_argument("--fmax", type=float, default=0.01)
    return p


def main():
    args = build_parser().parse_args()
    if args.mode == "calc":
        run_calc(args)
    elif args.mode == "plot":
        run_plot(args)
    else:
        run_calc(args)
        run_plot(args)


if __name__ == "__main__":
    main()
