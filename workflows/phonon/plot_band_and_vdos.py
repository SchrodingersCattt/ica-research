import numpy as np
import matplotlib.pyplot as plt
import yaml
import glob
import os

def parse_phonon_spectrum(data_file):
    with open(data_file) as f:
        lines = f.readlines()
        demarcation = [float(x) for x in lines[1].split('#')[-1].split()]

    q_segments, e_segments = [], []
    current_q, current_e = [], []

    for line in lines:
        stripped = line.strip()
        if not stripped:
            if current_q:
                q_segments.append(current_q)
                e_segments.append(current_e)
                current_q, current_e = [], []
            continue
        if stripped.startswith('#'):
            continue
        values = list(map(float, stripped.split()))
        current_q.append(values[0])
        current_e.append(values[1])

    if current_q:
        q_segments.append(current_q)
        e_segments.append(current_e)

    return q_segments, e_segments, demarcation

def regularize_bz(raw_lists):
    reg_bz = []
    for i in range(len(raw_lists) - 1):
        current_list = raw_lists[i]
        next_list = raw_lists[i + 1]
        if i == 0:
            reg_bz.append(current_list[0])
        if current_list[-1] == next_list[0]:
            reg_bz.append(current_list[-1])
        else:
            reg_bz.append(f"{current_list[-1]}|{next_list[0]}")
    if raw_lists:
        reg_bz.append(raw_lists[-1][-1])
    return reg_bz

def get_bz(file_path):
    with open(file_path, 'r') as file:
        data = yaml.safe_load(file)
    labels = data.get('labels', [])
    return regularize_bz(labels)

def plot_phonon_band(ax, _q, _e, ticks, labels, title, color):
    for q, e in zip(_q, _e):
        ax.plot(q, e, c=color, linewidth=0.5)
    ax.set_xlim(min(_q[0]), max(_q[-1]))
    ax.set_ylim(-5, 80)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, rotation=0)
    ax.axhline(y=0, color='k', linewidth=0.5, linestyle='--')
    for x in ticks:
        ax.axvline(x=x, color='k', linewidth=0.5)
    ax.set_ylabel("Frequency (THz)")
    # ax.set_title(title, fontsize=10)
    title = title.replace("Et4", r"Et$_4$")
    ax.text(0.05, 0.95, "  ", transform=ax.transAxes, fontsize=10)

def read_yaml(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        return yaml.safe_load(file)

def read_pdos_file(pdos_file):
    if not os.path.exists(pdos_file):
        raise FileNotFoundError(f"File not found: {pdos_file}")
    return np.loadtxt(pdos_file)

def extract_pdos_components(pdos_data, energy, elements, index_groups, color_dict):
    pdos_values, labels = [], []
    for element, indices in zip(elements, index_groups):
        label = format_label(element)
        pdos_sum = np.sum(pdos_data[:, indices], axis=1)
        pdos_values.append(pdos_sum)
        labels.append(label)
    return pdos_values, labels

def format_label(name):
    if 'N_dabco' in name:
        return r"N@H$_2$DABCO"
    elif 'N_ammonium' in name:
        return "N@ammonium"
    return name

def plot_pdos(ax, energy, pdos_vals, labels, color_dict):
    for label, data in zip(labels, pdos_vals):
        ax.plot(data, energy, label=label, c=color_dict.get(label, 'k'))
    ax.set_ylim(-5, 80)
    ax.set_xlim(left=0)
    ax.axhline(y=0, color='b', linewidth=0.5, linestyle='--')
    # add legend
    ax.legend(frameon=False, handlelength=1, fontsize=8, handletextpad=0.5, loc='upper right')
    ax.set_yticks([])
    ax.set_xticks([])
    # ax.text(0.05, 0.95, "PDOS", transform=ax.transAxes, fontsize=10)


def count_phonon_modes_in_band(energy, dos, w_min, w_max):
    mask = (energy >= w_min) & (energy <= w_max)
    dos_in_band = dos[mask]
    freq_in_band = energy[mask]
    if len(freq_in_band) < 2:
        return 0
    phonon_num = np.trapz(dos_in_band, freq_in_band)
    return phonon_num


# Main Plotting
if __name__ == '__main__':
    # System settings
    sys_dict = {
        'CA': ["Cu", "N"],
        'NEt4-CA': ['Cu', 'organic cation', 'azide'],
        'LA': ['Pb', 'N']
    }

    pdos_indices_dict = {
        'CA': [[0, 3, 1, 2], [11, 8, 7, 21, 6, 22, 17, 5, 4, 20, 25, 24, 12, 13, 19, 15, 16, 9, 10, 26, 27, 14, 18, 23]],
        'NEt4-CA': [
                [12, 10, 14, 5, 7, 8, 11, 1, 3, 13, 0, 2, 15, 4, 9, 6],
                [234, 201, 194, 166, 151, 221, 217, 143, 136, 146, 153, 156, 206, 138, 199, 198, 213, 175, 150, 179, 207, 126, 189, 237, 137, 209, 178, 220, 235, 239, 187, 157, 231, 152, 165, 142, 167, 212, 215, 174, 230, 186, 229, 225, 184, 144, 139, 159, 192, 238, 162, 200, 188, 180, 173, 147, 193, 129, 185, 216, 195, 127, 164, 128, 236, 202, 172, 208, 222, 223, 163, 145, 224, 228, 134, 135, 158, 214, 181, 203,
                
                124, 132, 204, 196, 125, 149, 232, 219, 160, 211, 182, 218, 154, 177, 148, 140, 190, 133, 170, 226, 183, 176, 210, 155, 161, 227, 205, 141, 197, 171, 233, 191,
                81, 82,  83, 84], # N@EtN
                [23, 67, 109, 59, 118, 29, 47, 40, 41, 42, 120, 73, 33, 70, 39, 51, 75, 121, 25, 72, 86, 28, 22, 99, 85, 17, 27, 46, 37, 68, 78, 77, 43, 168, 169, 100, 64, 16, 65, 49, 74, 52, 115, 98, 35, 119, 53, 96, 94, 80, 63, 36, 19, 76, 101, 48, 108, 89, 26, 112, 131, 21, 102, 20, 54, 92, 130, 31, 79, 87, 60, 88, 104, 18, 66, 71, 97, 32, 122, 93, 38, 69, 58, 24, 107, 62, 117, 50, 110, 30, 44, 55, 91, 45, 105, 103, 113, 56, 111,34, 57, 95, 114, 90, 106, 61, 123, 116] # N@azide
                ],
        'LA': [[0, 2, 6, 5, 4, 1, 9, 10, 11, 3, 8, 7], 
               [57, 56, 45, 17, 70, 27, 60, 22, 49, 65, 52, 51, 68, 67, 35, 48, 24, 83, 16, 14, 54, 26, 61, 30, 80, 62, 18, 36, 29, 41, 73, 69, 76, 53, 34, 25, 55, 74, 72, 78, 81, 58, 79, 20, 32, 75, 47, 19, 33, 43, 13, 82, 50, 38, 63, 40, 64, 15, 21, 23, 59, 42, 39, 66, 71, 31, 44, 28, 46, 37, 77, 12]]
    }

    ele_color_dict = {
        'H': '#cdcdcd',
        'C': '#848484',
        'organic cation': '#24B486',
        'N': '#8D76E2',
        'azide': '#8D76E2',
        "Cu": "#981212",
        "Pb": "#204E77"
    }
    sys_color_dict = {
        'CA': '#ff0000',
        'NEt4-CA': '#EA8A42',
        'LA': '#88C4EA'
    }

    band_data_files = glob.glob('output/*_band.dat')
    systems = [os.path.basename(f).split('_')[0] for f in band_data_files]


    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(8 * len(systems), 6))
    outer = GridSpec(1, len(systems), wspace=0.2)  # spacing between systems
    axes = []
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 12
    for i in range(len(systems)):
        inner = outer[i].subgridspec(1, 2, width_ratios=[3, 1], wspace=0.1)

        axes.append(plt.subplot(inner[0]))
        axes.append(plt.subplot(inner[1]))

    for i, system in enumerate(systems):
        band_file = f"output/{system}_band.dat"
        yaml_file = f"output/{system}_band.yaml"
        pdos_file = f"output/{system}_pdos.dat"

        bz_labels = get_bz(yaml_file)
        q, e, ticks = parse_phonon_spectrum(band_file)
        raw_indices = pdos_indices_dict[system]
        elements = sys_dict[system]
        pdos_data = read_pdos_file(pdos_file)
        energy = pdos_data[:, 0]
        # The groups are explicit zero-based atom indices, not contiguous
        # blocks. Converting them by group length mislabels most ICA azide N
        # atoms as organic and most C/H atoms as azide.
        parsed_indices = [[atom_index + 1 for atom_index in group] for group in raw_indices]
        pdos_vals, labels = extract_pdos_components(pdos_data, energy, elements, parsed_indices, ele_color_dict)

        ax_band = axes[2*i]
        ax_pdos = axes[2*i + 1]

        total_dos = np.sum(pdos_data[:, 1:], axis=1)
        w_at_gap = 10.0
        w_min = w_at_gap
        w_max = 2 * w_at_gap
        phonon_num = count_phonon_modes_in_band(energy, total_dos, w_min, w_max)
        print(f"{system} Phonon Modes w_c {w_min} 2w_c {w_max}: {phonon_num:.2f}")


        plot_phonon_band(ax_band, q, e, ticks, bz_labels, title=system, color=sys_color_dict.get(system, 'k'))
        plot_pdos(ax_pdos, energy, pdos_vals, labels, ele_color_dict)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.15)
    plt.savefig("band_and_dos.png", dpi=600)
    plt.savefig("Azides_phonon_and_dos.pdf")
