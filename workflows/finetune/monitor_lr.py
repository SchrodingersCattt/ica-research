import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import savgol_filter
import glob
import os
from tqdm import tqdm

STD_E = 1
STD_F = 1
STD_V = 1
lcurves = glob.glob('*/lcurve.out')


def read_file_skip_bad_lines(file_path):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            try:
                values = list(map(float, line.split()))
                data.append(values)
            except ValueError:
                continue
    return np.array(data)


def smooth(x, window_len):
    try:
        return savgol_filter(x, window_length=window_len, polyorder=1)
    except ValueError:
        return x


def get_eff(logfile=None):
    if logfile is None:
        return
    else:
        wall_times = []
        with open(logfile, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if "total wall time =" in line:
                time = float(line.split()[-2])
                wall_times.append(time)
        avg_time = np.mean(wall_times)
        return avg_time


plt.rcParams['font.size'] = 8
plt.rcParams['font.family'] = 'Arial'

n_curves = len(lcurves)
rows = int(np.ceil(n_curves / 1))
fig, axs = plt.subplots(1, rows, figsize=(3 * rows, 3))

# 检查 axs 是否为数组，如果不是则将其包装成数组
if rows == 1:
    axs = np.array([axs])

for idx, lcurvefile in enumerate(lcurves):
    try:
        lcurve = read_file_skip_bad_lines(lcurvefile)
        train_log = os.path.join(os.path.dirname(lcurvefile), 'train.log')
        if os.path.exists(train_log):
            eff = get_eff(train_log)
        else:
            eff = None

        if lcurve.shape[1] == 10:
            df = pd.DataFrame(lcurve, columns=['step', 'rmse_val', 'rmse_trn', 'rmse_e_val', 'rmse_e_trn',
                                               'rmse_f_val', 'rmse_f_trn', 'rmse_v_val', 'rmse_v_trn', 'lr'])
            has_val_columns = True
            has_v_columns = True
        elif lcurve.shape[1] == 8:
            df = pd.DataFrame(lcurve, columns=['step', 'rmse_val', 'rmse_trn', 'rmse_e_val', 'rmse_e_trn',
                                               'rmse_f_val', 'rmse_f_trn', 'lr'])
            has_val_columns = True
            has_v_columns = False
        elif lcurve.shape[1] == 6:
            df = pd.DataFrame(lcurve, columns=['step', 'rmse_trn', 'rmse_e_trn', 'rmse_f_trn', 'rmse_v_trn', 'lr'])
            has_val_columns = False
            has_v_columns = True
        else:
            print(f"Unexpected number of columns in {lcurvefile}: {lcurve.shape[1]}")
            continue

        # Applying the smoothing function to the RMSE values
        smoothed_rmse_e_trn = smooth(df['rmse_e_trn'], window_len=100)
        if has_val_columns:
            smoothed_rmse_e_val = smooth(df['rmse_e_val'], window_len=100)
        smoothed_rmse_f_trn = smooth(df['rmse_f_trn'], window_len=100)
        if has_val_columns:
            smoothed_rmse_f_val = smooth(df['rmse_f_val'], window_len=100)
        if has_v_columns:
            smoothed_rmse_v_trn = smooth(df['rmse_v_trn'], window_len=100)
            if has_val_columns:
                smoothed_rmse_v_val = smooth(df['rmse_v_val'], window_len=100)

    except Exception as e:
        print(lcurvefile)
        print(f"Error processing {lcurvefile}: {e}")
        continue

    ax1 = axs[idx]

    #ax1.axhline(STD_E, c='#cb4455', linestyle=':', label=r'std, $E$')
    ax1.plot(df['step'], smoothed_rmse_e_trn, c='#cb4455', linestyle='-', label=r'trn, $E$')
    if has_val_columns:
        ax1.plot(df['step'], smoothed_rmse_e_val, c='#cb4455', linestyle='-', label=r'val, $E$', alpha=0.5)
    #ax1.axhline(STD_F, c='#34ab56', linestyle=':', label=r'std, $F$')
    ax1.plot(df['step'], smoothed_rmse_f_trn, c='#34ab56', linestyle='-', label=r'trn, $F$')
    if has_val_columns:
        ax1.plot(df['step'], smoothed_rmse_f_val, c='#34ab56', linestyle='-', label=r'val, $F$', alpha=0.5)

    if has_v_columns:
        #ax1.axhline(STD_V, c='#456bfa', linestyle=':', label=r'std, $V$')
        ax1.plot(df['step'], smoothed_rmse_v_trn, c='#456bfa', linestyle='-', label=r'trn, $V$')
        if has_val_columns:
            ax1.plot(df['step'], smoothed_rmse_v_val, c='#456bfa', linestyle='-', label=r'val, $V$', alpha=0.5)

    ax1.text(
        0.05, 0.05,
        "\n".join(lcurvefile.split('/')[0].split('__')),
        transform=ax1.transAxes,
        horizontalalignment='left',
        verticalalignment='bottom',
        fontsize=8
    )
    if eff is not None:
        ax1.text(
            0.95, 0.05,
            f"Ava. wall time:\n{eff:.1f}s / 100steps",
            transform=ax1.transAxes,
            horizontalalignment='right',
            verticalalignment='bottom',
            fontsize=8
        )

    ax1.set_xlabel('Step')
    ax1.set_ylabel('Loss')

    ax1.set_yscale('log')
    ax1.set_xscale('log')

    ax1.set_xlim(1e2, 1e6)
    ax1.set_ylim(1e-3, 1e1)
    ax1.grid()
    if has_val_columns and has_v_columns:
        ncols = 3
    elif has_val_columns or has_v_columns:
        ncols = 3
    else:
        ncols = 1
    ax1.legend(fontsize=8, loc='upper left', ncols=ncols, handlelength=1.5, frameon=False, handletextpad=0.5, columnspacing=0.5)

plt.tight_layout()
plt.savefig('lcurve.png', dpi=300)