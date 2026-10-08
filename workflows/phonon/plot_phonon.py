import numpy as np
import matplotlib.pyplot as plt
import glob
import yaml
import os

def parse_phonon_spectrum(data_file):
    with open(data_file) as f:
        lines = f.readlines()
        demarcation = [float(x) for x in lines[1].split('#')[-1].split()]
    
    q_segments = []
    e_segments = []
    current_q = []
    current_e = []
    
    for line in lines:
        stripped = line.strip()
        if not stripped:  # Empty line: save current segment and reset
            if current_q:
                q_segments.append(current_q)
                e_segments.append(current_e)
                current_q = []
                current_e = []
            continue
        if stripped.startswith('#'):  # Skip comment lines
            continue
        values = list(map(float, stripped.split()))
        current_q.append(values[0])  # First column -> Q
        current_e.append(values[1])  # Second column -> E
    
    # Handle the last segment if no empty line at EOF
    if current_q:
        q_segments.append(current_q)
        e_segments.append(current_e)
    
    return q_segments, e_segments, demarcation

def regularize_bz(raw_lists):
    reg_bz = []
    for i in range(len(raw_lists) - 1):
        current_list = raw_lists[i]
        next_list = raw_lists[i + 1]
        # Add the first element of the current list
        if i == 0:
            reg_bz.append(current_list[0])
        # Merge the last element of current list with the first element of next list
        if current_list[-1] == next_list[0]:
            reg_bz.append(current_list[-1])
        else:
            merged_element = f"{current_list[-1]}|{next_list[0]}"
            reg_bz.append(merged_element)
    # Add the last element of the last list
    if raw_lists:
        reg_bz.append(raw_lists[-1][-1])
    return reg_bz

def get_bz(file_path):
    # *_band.yaml
    with open(file_path, 'r') as file:
        data = yaml.safe_load(file)
    labels = data.get('labels', [])
    return regularize_bz(labels)

def plot_phonon_spectrum(_q, _e, l, bz, **kwargs):
    color = kwargs.get('color', 'r')
    for q, e in zip(_q, _e):
        plt.plot(q, e, c=color, linewidth=0.5)
    plt.xlim(min(q), max(q))
    plt.ylim(-5, 80)
    plt.xticks(l, bz, rotation=0)
    # plt.title((band_data_file.split('/')[1].split('_')[0]).upper())
    for x in l:
        plt.axvline(x=x, color='k', linewidth=0.5) 
    plt.axhline(y=0, color='k', linewidth=0.5, linestyle='--')
    plt.xlabel("Brillouin Zone")
    plt.ylabel('Frequency (THz)')
    plt.tight_layout()

def save_for_origin(q_segments, e_segments, filename_prefix='phonon'):
    with open(f'{filename_prefix}.txt', 'w') as f:
        for q, e in zip(q_segments, e_segments):
            for x, y in zip(q, e):
                f.write(f"{x:.6f}\t{y:.6f}\n")
            f.write("\n")
    with open(f'{filename_prefix}_bz.txt', 'w') as f:
        f.write("Position\tLabel\n")
        for pos, label in zip(l, bz):
            f.write(f"{pos}\t{label}\n")


if __name__ == '__main__':
    band_data_files = glob.glob('output/*_band.dat')
    plt.figure(figsize=(12, 5))
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 12
    color_dict = {
        'CA': '#ff0000',
        'NEt4-CA': '#EA8A42',
        'LA': '#88C4EA'
    }
    # colors = ['#33afab', "#ab9900", '#a63399', '#3366a9']
    for i, band_data_file in enumerate(sorted(band_data_files), start=1):
        sys = band_data_file.split('/')[1].split('_')[0]
        color = color_dict.get(sys, 'k')
        bz = get_bz(band_data_file.replace('_band.dat', '_band.yaml'))
        _q, _e, l = parse_phonon_spectrum(band_data_file)
        save_for_origin(_q, _e, filename_prefix=band_data_file.split('/')[-1].replace('_band.dat',''))
        plt.subplot(1, 3, i)
        plot_phonon_spectrum(_q, _e, l, bz, color=color)
    plt.savefig('phonon.png', dpi=300)
