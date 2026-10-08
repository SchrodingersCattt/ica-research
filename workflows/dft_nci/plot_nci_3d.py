#!/usr/bin/env python3
"""
Nature-quality 3D NCI isosurface — PyVista v5.

Back to RDG=0.5 (correct NCI isovalue).
Focus: make the thin NCI patches visible by using proper rendering.
Atoms very light/wireframe, surface bold and colorful.
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree
from skimage.measure import marching_cubes
import pyvista as pv
from matplotlib.colors import ListedColormap
import os, time

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
COLOR_FILE = os.path.join(SCRIPT_DIR, "func1.cub")
ISO_FILE   = os.path.join(SCRIPT_DIR, "func2.cub")

ISOVALUE = 0.5              # Standard NCI isovalue
MESH_SMOOTH = 50            # Light smoothing only

CMIN = -0.035
CMAX =  0.02

BOHR_TO_ANG = 0.529177

ELEMENT_DATA = {
    1:  ("H",  [0.90, 0.90, 0.90]),
    6:  ("C",  [0.45, 0.45, 0.45]),
    7:  ("N",  [0.35, 0.50, 0.92]),
    8:  ("O",  [0.90, 0.15, 0.15]),
    29: ("Cu", [0.85, 0.60, 0.30]),
}

ATOM_RADIUS = {'Cu': 0.30, 'N': 0.12, 'C': 0.12, 'H': 0.06}
BOND_RADIUS = 0.05
BOND_CUTOFF = 2.2
BOND_CUTOFF_METAL = 2.8

OUTPUT_PNG = os.path.join(SCRIPT_DIR, "NCI_3d_isosurface.png")
FIG_W, FIG_H = 2000, 2800


def read_cube(filepath):
    print(f"  Reading {os.path.basename(filepath)} ...")
    t0 = time.time()
    with open(filepath, 'r') as f:
        f.readline(); f.readline()
        parts = f.readline().split()
        natoms = int(parts[0])
        origin = np.array([float(x) for x in parts[1:4]])
        npts = np.zeros(3, dtype=int)
        axes = np.zeros((3, 3))
        for i in range(3):
            p = f.readline().split()
            npts[i] = int(p[0])
            axes[i] = [float(p[1]), float(p[2]), float(p[3])]
        atoms = []
        for _ in range(natoms):
            p = f.readline().split()
            atoms.append((int(p[0]), float(p[1]), float(p[2]), float(p[3]), float(p[4])))
        vals = []
        for line in f:
            vals.extend(line.split())
    data = np.array(vals, dtype=np.float64).reshape(npts[0], npts[1], npts[2])
    print(f"    {npts[0]}x{npts[1]}x{npts[2]}, [{data.min():.3f}, {data.max():.3f}] ({time.time()-t0:.1f}s)")
    return atoms, origin, axes, npts, data


def find_bonds(xyz, Z_arr):
    tree = cKDTree(xyz)
    bonds = set()
    for i, j in tree.query_pairs(r=BOND_CUTOFF):
        zi, zj = Z_arr[i], Z_arr[j]
        if zi == 1 and zj == 1: continue
        if (zi == 1 or zj == 1) and np.linalg.norm(xyz[i]-xyz[j]) > 1.3: continue
        if zi == 29 and zj == 29: continue
        bonds.add((min(i,j), max(i,j)))
    for i, j in tree.query_pairs(r=BOND_CUTOFF_METAL):
        zi, zj = Z_arr[i], Z_arr[j]
        if {zi, zj} == {29, 7}:
            bonds.add((min(i,j), max(i,j)))
    return list(bonds)


def make_nci_cmap(n=256):
    cpts = [
        (0.00, (0,   20,  230)),
        (0.20, (0,   120, 250)),
        (0.40, (0,   220, 100)),
        (0.55, (30,  230, 30)),
        (0.70, (240, 220, 0)),
        (0.85, (255, 100, 0)),
        (1.00, (230, 10,  10)),
    ]
    colors = np.zeros((n, 3))
    for k in range(n):
        t = k / (n - 1)
        for idx in range(len(cpts) - 1):
            t0, c0 = cpts[idx]; t1, c1 = cpts[idx + 1]
            if t0 <= t <= t1:
                f = (t - t0) / (t1 - t0) if t1 > t0 else 0
                colors[k] = [c0[i] + f * (c1[i] - c0[i]) for i in range(3)]
                break
    return ListedColormap(colors / 255.0)


def main():
    print("=" * 60)
    print("NCI 3D v5 — RDG=0.5, clean rendering")
    print("=" * 60)

    atoms_raw, orig, axes, npts, color_data = read_cube(COLOR_FILE)
    _, _, _, _, iso_data = read_cube(ISO_FILE)

    x = orig[0] + np.arange(npts[0]) * axes[0, 0]
    y = orig[1] + np.arange(npts[1]) * axes[1, 1]
    z = orig[2] + np.arange(npts[2]) * axes[2, 2]
    xa, ya, za = x * BOHR_TO_ANG, y * BOHR_TO_ANG, z * BOHR_TO_ANG
    sp = (abs(axes[0,0])*BOHR_TO_ANG, abs(axes[1,1])*BOHR_TO_ANG, abs(axes[2,2])*BOHR_TO_ANG)

    # Marching cubes
    print(f"\n  Marching cubes at RDG={ISOVALUE} ...")
    verts, faces, normals, _ = marching_cubes(iso_data, level=ISOVALUE, spacing=sp)
    verts[:, 0] += xa[0]; verts[:, 1] += ya[0]; verts[:, 2] += za[0]
    print(f"    {len(verts):,} verts, {len(faces):,} faces")

    # Build mesh with normals for proper two-sided rendering
    faces_pv = np.column_stack([np.full(len(faces), 3), faces]).ravel()
    mesh = pv.PolyData(verts, faces_pv)

    # Color mapping
    print("  Color mapping ...")
    interp = RegularGridInterpolator(
        (xa, ya, za), color_data, method='linear', bounds_error=False, fill_value=0.0
    )
    mesh["sign_l2_rho"] = np.clip(interp(np.array(mesh.points)), CMIN, CMAX)

    # Light smoothing
    if MESH_SMOOTH > 0:
        print(f"  Taubin smoothing ({MESH_SMOOTH} iters) ...")
        mesh = mesh.smooth_taubin(n_iter=MESH_SMOOTH, pass_band=0.1)
        mesh["sign_l2_rho"] = np.clip(interp(np.array(mesh.points)), CMIN, CMAX)

    # Compute normals for both sides
    mesh = mesh.compute_normals(auto_orient_normals=True)

    # Atoms
    atom_xyz, atom_Z = [], []
    for Z, ch, ax, ay, az in atoms_raw:
        if Z == 1: continue  # no H
        atom_xyz.append([ax * BOHR_TO_ANG, ay * BOHR_TO_ANG, az * BOHR_TO_ANG])
        atom_Z.append(Z)
    atom_xyz = np.array(atom_xyz)
    atom_Z = np.array(atom_Z)

    # ═══ Render ═══
    print("\n  Rendering ...")
    pv.OFF_SCREEN = True
    p = pv.Plotter(off_screen=True, window_size=[FIG_W, FIG_H])
    p.set_background([0.98, 0.98, 0.98])  # very light gray bg

    nci_cmap = make_nci_cmap(256)

    # NCI surface — bold, fully opaque, two-sided
    p.add_mesh(
        mesh,
        scalars="sign_l2_rho",
        cmap=nci_cmap,
        clim=[CMIN, CMAX],
        opacity=1.0,           # FULLY OPAQUE
        smooth_shading=True,
        show_edges=False,
        backface_culling=False,  # show both sides
        specular=0.5,
        specular_power=40,
        diffuse=0.7,
        ambient=0.3,
        show_scalar_bar=True,
        scalar_bar_args=dict(
            title="sign(λ₂)ρ (a.u.)",
            title_font_size=20,
            label_font_size=15,
            n_labels=5,
            position_x=0.82,
            position_y=0.18,
            width=0.06,
            height=0.55,
            fmt="%.3f",
            font_family="arial",
            color='black',
        ),
    )

    # Atoms — very faint, wireframe-like
    print("  Adding atoms ...")
    for Z_val in np.unique(atom_Z):
        mask = atom_Z == Z_val
        pts = atom_xyz[mask]
        sym, col = ELEMENT_DATA.get(Z_val, ("?", [0.5, 0.5, 0.5]))[:2]
        r = ATOM_RADIUS.get(sym, 0.10)
        # Lighter colors
        col_faint = [min(1.0, c * 0.5 + 0.5) for c in col]
        
        cloud = pv.PolyData(pts)
        sphere = pv.Sphere(radius=r, theta_resolution=14, phi_resolution=14)
        glyphs = cloud.glyph(geom=sphere, orient=False, scale=False)
        p.add_mesh(
            glyphs, color=col_faint,
            smooth_shading=True,
            opacity=0.35,
            specular=0.1, diffuse=0.4, ambient=0.5,
        )

    # Bonds — very faint
    print("  Adding bonds ...")
    bonds = find_bonds(atom_xyz, atom_Z)
    print(f"    {len(bonds)} bonds")
    all_pts, all_lines = [], []
    idx = 0
    for i, j in bonds:
        all_pts.extend([atom_xyz[i], atom_xyz[j]])
        all_lines.append([2, idx, idx + 1])
        idx += 2
    if all_pts:
        lines = pv.PolyData(np.array(all_pts), lines=np.hstack(all_lines))
        tubes = lines.tube(radius=BOND_RADIUS, n_sides=10)
        p.add_mesh(
            tubes, color=[0.75, 0.75, 0.75],
            smooth_shading=True,
            opacity=0.30,
            specular=0.05,
        )

    # Camera
    cx = (xa[0] + xa[-1]) / 2
    cy = (ya[0] + ya[-1]) / 2
    cz = (za[0] + za[-1]) / 2

    # Front view along x, y up
    p.camera_position = [
        (cx + 80, cy, cz + 6),
        (cx, cy, cz),
        (0, 1, 0),
    ]
    p.camera.parallel_projection = True
    p.camera.parallel_scale = (ya[-1] - ya[0]) * 0.55

    # Lighting
    p.remove_all_lights()
    p.add_light(pv.Light(
        position=(cx + 50, cy + 40, cz + 70),
        focal_point=(cx, cy, cz),
        intensity=1.0,
        light_type='scenelight',
    ))
    p.add_light(pv.Light(
        position=(cx - 30, cy - 20, cz + 30),
        focal_point=(cx, cy, cz),
        intensity=0.4,
        light_type='scenelight',
    ))

    print(f"\n  Saving -> {OUTPUT_PNG}")
    p.screenshot(OUTPUT_PNG, transparent_background=False)
    p.close()

    print("Done!")


if __name__ == "__main__":
    main()
