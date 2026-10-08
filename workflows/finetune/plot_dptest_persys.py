import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from time import time
from scipy.stats import gaussian_kde

np.random.seed(42)
plt.style.use('seaborn-v0_8-paper') 
plt.rcParams['font.family'] = "Arial"

# Define color map anchors
# color_anchor = ["#44BB89", "#8844ff", "#EE4455"]
color_anchor = ["#CCCCCC", "#8484cc", "#0033F8"]
cmap = LinearSegmentedColormap.from_list('colorbar', color_anchor)

plt.rcParams['xtick.major.pad'] = 0.5  # Padding for x-axis major ticks
plt.rcParams['ytick.major.pad'] = 0.5  # Padding for y-axis major ticks

def calc_rmse(gt, pred):
    if isinstance(gt, list):
        gt = np.array(gt)
        pred = np.array(pred)
    return np.sqrt(np.mean((gt - pred)**2))


def filter_data(data1, data2, num_points):
    """Randomly select a specified number of points from the data."""
    if (num_points == 1) or (not num_points):
        return data1, data2
    else:
        indices = np.random.choice(len(data1), size=num_points, replace=False)
        return data1[indices], data2[indices]

def parse_e(e_file, mode='direct'):
    if mode == 'direct':
        e = np.loadtxt(e_file)
        e_dft = e[:, 0]
        e_dp = e[:, 1]
        return e_dft, e_dp
    elif mode == 'eliminate_lowest':
        with open(e_file, 'r') as file:
            content = file.read()
            datasets = content.strip().split('#')[1:]

        e_dft = []
        e_dp = []
        for data_blob in datasets:
            dataset_name = data_blob.split(':', 1)[0].strip()
            print(f"Processing dataset: {dataset_name}")

            lines = data_blob.split('\n')[1:]  # Skip the line with dataset name
            energy_dft = []
            energy_dp = []

            for line in lines:
                if not line.strip(): continue
                if line.startswith("#"): break
                values = line.split()
                energy_dft.append(float(values[0]))  # DFT energies
                energy_dp.append(float(values[1]))    # DP energies
            
            energy_dft = np.array([ee - min(energy_dft) for ee in energy_dft])
            energy_dp = np.array([ee - min(energy_dp) for ee in energy_dp])

            e_dft.extend(energy_dft)
            e_dp.extend(energy_dp)

        e_dft = np.array(e_dft)  
        e_dp = np.array(e_dp)     
        return e_dft, e_dp
    elif mode == 'eliminate bias':
        pass

def load_and_process_data(e_file, f_file, v_file=None, parse_e_mode=None):
    start = time()
    e = np.loadtxt(e_file)
    f = np.loadtxt(f_file)
    e_dft, e_dp = parse_e(e_file, mode=parse_e_mode)

    ff_dft = np.stack((f[:, 0], f[:, 1], f[:, 2]), axis=0).flatten()
    ff_dp = np.stack((f[:, 3], f[:, 4], f[:, 5]), axis=0).flatten()
    if v_file is not None:
        v = np.loadtxt(v_file)
        vv_dft = np.stack((v[:, 0], v[:, 1], v[:, 2],
                           v[:, 3], v[:, 4], v[:, 5],
                           v[:, 6], v[:, 7], v[:, 8]), axis=0).flatten()
        vv_dp = np.stack((v[:, 9], v[:, 10], v[:, 11],
                          v[:, 12], v[:, 13], v[:, 14],
                          v[:, 15], v[:, 16], v[:, 17],), axis=0).flatten()
        end = time()
        print(end - start)
        return e_dft, e_dp, ff_dft, ff_dp, vv_dft, vv_dp
    else:
        print('No v file input')
        end = time()
        print(end - start)
        return e_dft, e_dp, ff_dft, ff_dp


def plot_comparison(e_file, f_file, v_file, parse_e_mode):
    e_dft, e_dp, f_data_dft, f_data_dp, v_data_dft, v_data_dp = load_and_process_data(e_file, f_file, v_file, parse_e_mode)
    if not parse_e_mode == 'direct':
        e_dft, e_dp, _, _, _, _ = load_and_process_data(e_file, f_file, v_file, 'direct')
    start = time()
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(5, 2))
    plt.rcParams['font.size'] = 12
    plt.rcParams['font.family'] = 'Arial'
    rmse_e, rmse_f, rmse_v = calc_rmse(e_dft, e_dp), calc_rmse(f_data_dft, f_data_dp), calc_rmse(v_data_dft, v_data_dp)

    print("rmse e", calc_rmse(e_dft, e_dp))
    print("rmse f", calc_rmse(f_data_dft, f_data_dp))
    print("rmse v", calc_rmse(v_data_dft, v_data_dp))

    # Filter the data if necessary
    max_points = 5000

    # Energy plot
    if len(e_dft) > max_points:
        e_dft, e_dp = filter_data(e_dft, e_dp, max_points)
    min_e = min(e_dft.min(), e_dp.min()) * 1.2
    max_e = max(e_dft.max(), e_dp.max()) * 0.8

    # Calculate density for coloring
    xy = np.vstack([e_dft, e_dp])
    density = gaussian_kde(xy)(xy)

    scatter1 = ax1.scatter(e_dft, e_dp, s=5, c=density, cmap=cmap)
    ax1.plot([min_e, max_e], [min_e, max_e], 'k--', linewidth=0.5,  label='DP=DFT')
    ax1.set_xlim(min_e, max_e)
    ax1.set_ylim(min_e, max_e)
    ax1.set_xlabel(r'$E_\mathrm{DFT} \ (\mathrm{eV/atom})$', labelpad=0.2)
    ax1.set_ylabel(r'$E_\mathrm{DP} \ (\mathrm{eV/atom})$', labelpad=0.2)
    #ax1.grid()
    ax1.text(0.05, 0.8, f"RMSE:\n{rmse_e*1000:.1f} meV/atom", transform=ax1.transAxes, fontsize=8)
    ax1.legend(handlelength=1, frameon=False, loc='lower right')
    cbar1 = plt.colorbar(scatter1, ax=ax1, pad=0.0)
    cbar1.ax.tick_params(length=0)  # Remove tick marks
    cbar1.set_ticks([])             # Remove tick labels

    # Force plot
    if len(f_data_dft) > max_points:
        f_data_dft, f_data_dp = filter_data(f_data_dft, f_data_dp, max_points)

    min_f = min(f_data_dft.min(), f_data_dp.min()) * 1.2
    max_f = max(f_data_dft.max(), f_data_dp.max()) * 1.2

    xy_f = np.vstack([f_data_dft, f_data_dp])
    density_f = gaussian_kde(xy_f)(xy_f)

    scatter2 = ax2.scatter(f_data_dft, f_data_dp, s=5, c=density_f, cmap=cmap)
    ax2.plot([min_f, max_f], [min_f, max_f], 'k--', linewidth=0.5,  label='DP=DFT')
    ax2.set_xlabel(r'$F_\mathrm{DFT} \ (\mathrm{eV/\AA})$', labelpad=0.2)
    ax2.set_ylabel(r'$F_\mathrm{DP} \ (\mathrm{eV/\AA})$', labelpad=0.2)
    ax2.set_xlim(min_f, max_f)
    ax2.set_ylim(min_f, max_f)
    #ax2.grid()
    ax2.text(0.05, 0.8, f"RMSE:\n{rmse_f*1000:.1f} meV/"+r"$\mathrm{\AA}$", transform=ax2.transAxes, fontsize=8)
    ax2.legend(handlelength=1, frameon=False, loc='lower right')
    cbar2 = plt.colorbar(scatter2, ax=ax2, pad=0.0)
    cbar2.ax.tick_params(length=0)  # Remove tick marks
    cbar2.set_ticks([])             # Remove tick labels

    # Virial plot
    if len(v_data_dft) > max_points:
        v_data_dft, v_data_dp = filter_data(v_data_dft, v_data_dp, max_points)

    min_v = min(v_data_dft.min(), v_data_dp.min()) * 1.2
    max_v = max(v_data_dft.max(), v_data_dp.max()) * 1.2

    xy_v = np.vstack([v_data_dft, v_data_dp])
    print(xy_v)
    density_v = gaussian_kde(xy_v)(xy_v)

    scatter3 = ax3.scatter(v_data_dft, v_data_dp, s=5, c=density_v, cmap=cmap)
    ax3.plot([min_v, max_v], [min_v, max_v], 'k--', linewidth=0.5,  label='DP=DFT')
    ax3.set_xlabel(r'$V_\mathrm{DFT} \ (\mathrm{eV/atom})$', labelpad=0.2)
    ax3.set_ylabel(r'$V_\mathrm{DP} \ (\mathrm{eV/atom})$', labelpad=0.2)
    ax3.set_xlim(min_v, max_v)
    ax3.set_ylim(min_v, max_v)
    #ax3.grid()
    ax3.text(0.05, 0.8, f"RMSE:\n{rmse_v*1000:.1f} meV/atom", transform=ax3.transAxes, fontsize=8)
    ax3.legend(handlelength=1, frameon=False, loc='lower right')
    cbar3 = plt.colorbar(scatter3, ax=ax3, pad=0.0)
    cbar3.ax.tick_params(length=0)  # Remove tick marks
    cbar3.set_ticks([])             # Remove tick labels

    plt.tight_layout(pad=0.8)
    plt.savefig(f"efv_comparison.png", dpi=300)
    plt.savefig(f"efv_comparison.eps")
    plt.close()

def plot_hist(e_file, f_file, v_file, parse_e_mode=None):
    e_dft, e_dp, f_data_dft, f_data_dp, v_data_dft, v_data_dp = load_and_process_data(e_file, f_file, v_file, parse_e_mode)
    plt.rcParams['font.size'] = 12
    plt.figure(figsize=(8, 4))
    plt.subplot(131)
    plt.hist(e_dft, bins=100, alpha=0.5, label='DFT')
    plt.hist(e_dp, bins=100, alpha=0.5, label='DP')
    plt.xlabel(r'$E_\mathrm{DFT} \ (\mathrm{eV/atom})$', labelpad=0.2)
    plt.ylabel('Counts')
    plt.legend(frameon=False)
    plt.subplot(132)
    plt.hist(f_data_dft, bins=100, alpha=0.5, label='DFT')
    plt.hist(f_data_dp, bins=100, alpha=0.5, label='DP')
    plt.xlabel(r'$F_\mathrm{DFT} \ (\mathrm{eV/\AA})$', labelpad=0.2)
    plt.legend(frameon=False)
    plt.subplot(133)
    plt.hist(v_data_dft, bins=100, alpha=0.5, label='DFT')
    plt.hist(v_data_dp, bins=100, alpha=0.5, label='DP')
    plt.xlabel(r'$V_\mathrm{DFT} \ (\mathrm{eV/atom})$', labelpad=0.2)
    plt.legend(frameon=False)
    plt.tight_layout()

    plt.savefig("efv_hist.png", dpi=300)

def main():
    e_file = 'test.e_peratom.out'
    f_file = 'test.f.out'
    v_file = 'test.v_peratom.out'
    
    plot_comparison(e_file=e_file, f_file=f_file, v_file=v_file, parse_e_mode='direct')
    plot_hist(e_file=e_file, f_file=f_file, v_file=v_file, parse_e_mode='eliminate_lowest')
if __name__ == "__main__":
    main()
