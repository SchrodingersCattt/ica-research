import os
import yaml
import numpy as np
import matplotlib.pyplot as plt


def read_yaml(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        return yaml.safe_load(file)


def read_pdos_file(pdos_file):
    if not os.path.exists(pdos_file):
        raise FileNotFoundError(f"File not found: {pdos_file}")
    return np.loadtxt(pdos_file)


def plot_pdos_spectrum(e, pdos, labels, colors, title, **kwargs):
    for label, data in zip(labels, pdos):
        plt.plot(
            e, data, 
            label="        ", 
            c=colors[label])

    plt.plot([max(e), 100], [0, 0], c="k", linewidth=0.5)
    plt.xlim(0, 80)
    plt.xlabel(r'Frequency (THz)')
    plt.ylabel('Phonon DOS')
    plt.grid(True)
    # plt.title(title, fontsize=10)
    plt.legend(ncols=3, handlelength=1.5, frameon=False)
    plt.tight_layout()


def extract_pdos_components(pdos_data, energy, elements, index_groups, color_dict):
    pdos_values = []
    labels = []

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


def save_pdos_for_origin(filename, energy, pdos_list, labels):
    header = "Frequency(THz)\t" + "\t".join(labels)
    data = np.column_stack([energy] + pdos_list)
    np.savetxt(filename, data, fmt="%.6f", delimiter='\t', header=header, comments='')

def main():
    sys_dict = {
        'CA': ["Cu", "azide"],
        'NEt4-CA': [
            'Cu', 
            'organic cation', 
            'azide'],
        'LA': ['Pb', 'azide']
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

    color_dict = {
        'H': '#cdcdcd',
        'C': '#848484',
        'organic cation': '#24B486',
        'N': '#8D76E2',
        'azide': '#8D76E2',
        "Cu": "#981212",
        "Pb": "#204E77"
    }

    plt.figure(figsize=(6, 6))
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 12

    for idx, (sys, elements) in enumerate(sorted(sys_dict.items()), start=1):
        pdos_file = f"output/{sys}_pdos.dat"
        raw_indices = pdos_indices_dict[sys]
        parsed_indices = [[i + 1 for i in group] for group in raw_indices]

        pdos_data = read_pdos_file(pdos_file)
        energy = pdos_data[:, 0]
        pdos_vals, labels = extract_pdos_components(pdos_data, energy, elements, parsed_indices, color_dict)

        # Pad data with zeros up to 80 THz if needed
        if energy[-1] < 80:
            # Create new energy array from 0 to 80 with same step as original data
            energy_step = energy[1] - energy[0]
            new_energy_end = int(80/energy_step) * energy_step
            new_energy_points = np.arange(energy[-1] + energy_step, new_energy_end + energy_step, energy_step)
            
            # Extend energy array
            energy = np.concatenate([energy, new_energy_points])
            
            # Extend each pdos component with zeros
            new_pdos_vals = []
            for pdos_component in pdos_vals:
                zeros_padding = np.zeros(len(new_energy_points))
                new_pdos_component = np.concatenate([pdos_component, zeros_padding])
                new_pdos_vals.append(new_pdos_component)
            pdos_vals = new_pdos_vals


        plt.subplot(3, 1, idx)
        plot_pdos_spectrum(energy, pdos_vals, labels, color_dict, title=sys)      
        save_pdos_for_origin(f"output/{sys}_pdos.txt", energy, pdos_vals, labels)


    plt.savefig("vib_pdos.png", dpi=300)


if __name__ == '__main__':
    main()