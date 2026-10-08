#!/usr/bin/env python3
"""
Nature-quality NCI/IRI 2D scatter plot: RDG vs sign(λ₂)ρ
Uses 2D histogram with density-based coloring for publication quality.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import gaussian_kde
import os

# ──────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────
DATA_FILE = os.path.join(os.path.dirname(__file__), "output.txt")
OUTPUT_PNG = os.path.join(os.path.dirname(__file__), "NCI_scatter_nature.png")
OUTPUT_PDF = os.path.join(os.path.dirname(__file__), "NCI_scatter_nature.pdf")

# Plot range
X_RANGE = (-0.05, 0.05)   # sign(λ₂)ρ range (a.u.)
Y_RANGE = (0.0, 2.0)      # RDG / IRI range
NBINS = 400               # histogram resolution

# ──────────────────────────────────────────────
# Load data (columns: x, y, z, sign(λ₂)ρ, RDG)
# ──────────────────────────────────────────────
print("Loading scatter data...")
data = np.loadtxt(DATA_FILE)
sign_lambda2_rho = data[:, 3]  # 4th column
rdg = data[:, 4]               # 5th column

# Filter to plot range
mask = (
    (sign_lambda2_rho >= X_RANGE[0]) & (sign_lambda2_rho <= X_RANGE[1]) &
    (rdg >= Y_RANGE[0]) & (rdg <= Y_RANGE[1])
)
x = sign_lambda2_rho[mask]
y = rdg[mask]
print(f"  Total points: {len(sign_lambda2_rho)}, in range: {len(x)}")

# ──────────────────────────────────────────────
# Create 2D histogram for density coloring
# ──────────────────────────────────────────────
print("Computing 2D density histogram...")
H, xedges, yedges = np.histogram2d(x, y, bins=NBINS, range=[X_RANGE, Y_RANGE])
H = H.T  # transpose for imshow orientation

# Log-scale density for better contrast
H_log = np.log10(H + 1)

# ──────────────────────────────────────────────
# Nature-style BWG colormap (Blue-White-Green-Red)
# Matches the standard NCI coloring convention:
#   Blue (negative) = attractive / H-bond
#   Green (near zero) = vdW
#   Red (positive) = repulsive / steric
# ──────────────────────────────────────────────
# Custom colormap: white background -> density gradient
density_cmap = LinearSegmentedColormap.from_list(
    'nci_density',
    [
        (1.0, 1.0, 1.0),     # white (background)
        (0.7, 0.85, 1.0),    # light blue
        (0.2, 0.4, 0.8),     # blue
        (0.1, 0.1, 0.6),     # dark blue
        (0.05, 0.05, 0.3),   # very dark blue
        (0.0, 0.0, 0.0),     # black (highest density)
    ]
)

# ──────────────────────────────────────────────
# Alternative: color each point by its sign(λ₂)ρ value (NCI convention)
# ──────────────────────────────────────────────
# Standard NCI colormap: Blue -> Green -> Red
nci_cmap = LinearSegmentedColormap.from_list(
    'nci_bgr',
    [
        (0.0, 0.0, 1.0),     # blue  (attractive, negative sign(λ₂)ρ)
        (0.0, 0.8, 0.0),     # green (vdW, near zero)
        (1.0, 0.0, 0.0),     # red   (repulsive, positive sign(λ₂)ρ)
    ]
)

# ──────────────────────────────────────────────
# Figure setup — Nature style
# ──────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.2,
    'ytick.major.width': 1.2,
    'xtick.major.size': 5,
    'ytick.major.size': 5,
    'xtick.minor.size': 3,
    'ytick.minor.size': 3,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
})

# ──────────────────────────────────────────────
# Plot 1: Density-colored scatter (colored by point density)
# ──────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(4.5, 3.5), dpi=300)

# Use imshow for the 2D histogram
extent = [X_RANGE[0], X_RANGE[1], Y_RANGE[0], Y_RANGE[1]]
im = ax.imshow(
    H_log, origin='lower', extent=extent,
    aspect='auto', cmap=density_cmap,
    interpolation='gaussian',
    vmin=0, vmax=H_log.max()
)

ax.set_xlabel(r'sign($\lambda_2$)$\rho$ (a.u.)', fontsize=12, labelpad=4)
ax.set_ylabel('Reduced density gradient', fontsize=12, labelpad=4)
ax.set_xlim(X_RANGE)
ax.set_ylim(Y_RANGE)

# Fine tick settings
ax.xaxis.set_major_locator(plt.MultipleLocator(0.01))
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.005))
ax.yaxis.set_major_locator(plt.MultipleLocator(0.5))
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.1))

# Colorbar
cbar = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
cbar.set_label('log$_{10}$(count + 1)', fontsize=10)
cbar.ax.tick_params(labelsize=8)

fig.tight_layout(pad=0.5)
fig.savefig(OUTPUT_PNG, dpi=600, bbox_inches='tight', transparent=False)
fig.savefig(OUTPUT_PDF, dpi=600, bbox_inches='tight', transparent=False)
print(f"  Saved: {OUTPUT_PNG}")
print(f"  Saved: {OUTPUT_PDF}")
plt.close(fig)

# ──────────────────────────────────────────────
# Plot 2: NCI-colored scatter (colored by sign(λ₂)ρ value — the standard NCI plot)
# ──────────────────────────────────────────────
OUTPUT_NCI_PNG = os.path.join(os.path.dirname(__file__), "NCI_scatter_colored.png")
OUTPUT_NCI_PDF = os.path.join(os.path.dirname(__file__), "NCI_scatter_colored.pdf")

fig2, ax2 = plt.subplots(figsize=(4.5, 3.5), dpi=300)

# Subsample for reasonable scatter performance
N_SAMPLE = min(200000, len(x))
rng = np.random.default_rng(42)
idx = rng.choice(len(x), size=N_SAMPLE, replace=False)
xs, ys = x[idx], y[idx]

# Sort by density (plot high-density points on top)
sc = ax2.scatter(
    xs, ys, c=xs, cmap=nci_cmap,
    s=0.3, alpha=0.6, edgecolors='none',
    vmin=X_RANGE[0], vmax=X_RANGE[1],
    rasterized=True
)

ax2.set_xlabel(r'sign($\lambda_2$)$\rho$ (a.u.)', fontsize=12, labelpad=4)
ax2.set_ylabel('Reduced density gradient', fontsize=12, labelpad=4)
ax2.set_xlim(X_RANGE)
ax2.set_ylim(Y_RANGE)

ax2.xaxis.set_major_locator(plt.MultipleLocator(0.01))
ax2.xaxis.set_minor_locator(plt.MultipleLocator(0.005))
ax2.yaxis.set_major_locator(plt.MultipleLocator(0.5))
ax2.yaxis.set_minor_locator(plt.MultipleLocator(0.1))

# NCI colorbar with annotations
cbar2 = fig2.colorbar(sc, ax=ax2, shrink=0.85, pad=0.02)
cbar2.set_label(r'sign($\lambda_2$)$\rho$ (a.u.)', fontsize=10)
cbar2.ax.tick_params(labelsize=8)

# Add text annotations on colorbar
cbar2.ax.text(0.5, 0.02, 'Attractive', transform=cbar2.ax.transAxes,
              ha='center', va='bottom', fontsize=7, color='blue', fontweight='bold')
cbar2.ax.text(0.5, 0.50, 'vdW', transform=cbar2.ax.transAxes,
              ha='center', va='center', fontsize=7, color='green', fontweight='bold')
cbar2.ax.text(0.5, 0.98, 'Repulsive', transform=cbar2.ax.transAxes,
              ha='center', va='top', fontsize=7, color='red', fontweight='bold')

fig2.tight_layout(pad=0.5)
fig2.savefig(OUTPUT_NCI_PNG, dpi=600, bbox_inches='tight', transparent=False)
fig2.savefig(OUTPUT_NCI_PDF, dpi=600, bbox_inches='tight', transparent=False)
print(f"  Saved: {OUTPUT_NCI_PNG}")
print(f"  Saved: {OUTPUT_NCI_PDF}")
plt.close(fig2)

print("Done! 2D scatter plots generated.")
