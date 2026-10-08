import os
import glob
import dpdata
import numpy as np
import pandas as pd
from deepmd.infer.deep_eval import DeepEval
from pathlib import Path
from tqdm import tqdm
from typing import List, Optional, Tuple
import warnings
warnings.filterwarnings("ignore")


all_type_map = [
    "H", "He", 
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", 
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", 
    "Ga", "Ge", "As", "Se", "Br", "Kr", 

    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", 
    "In", "Sn", "Sb", "Te", "I", "Xe",

    "Cs", "Ba", 
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", 
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", 
    "Tl", "Pb", "Bi", "Po", "At", "Rn",

    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
]


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

def get_step_prog(lcurve_file=None):
    if lcurve_file is None:
        return
    else:
        lcurve = np.loadtxt(lcurve_file)
        steps = lcurve[:, 0]
        return steps[-1]

def rmse(x, y):
    return np.sqrt(np.mean((x - y) ** 2))

class DPPTPredict:
    def load_model(self, model: Path):
        self.dp = DeepEval(model)
        self.type_map = self.dp.get_type_map()

    def evaluate(self,
                 coord: np.ndarray,
                 cell: Optional[np.ndarray],
                 atype: List[int]
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        coord = coord.reshape([1, -1, 3])
        if cell is not None:
            cell = cell.reshape([1, 3, 3])
        atype = atype.reshape([1, -1])
        e, f, v = self.dp.eval(coord, cell, atype, infer_batch_size=1)
        return e.reshape([1])[0], f.reshape([-1, 3]), v.reshape([3, 3])

    def monitor_rmse(self, input_path):
        infer_energies, infer_forces, infer_virials = [], [], []
        gt_energies, gt_forces, gt_virials = [], [], []
        for f in Path(input_path).rglob("type.raw"):
            sys = f.parent
            d = dpdata.MultiSystems()
            mixed_type = len(list(sys.glob("*/real_atom_types.npy"))) > 0
            if mixed_type:
                d.load_systems_from_file(sys, fmt="deepmd/npy/mixed")
            else:
                k = dpdata.LabeledSystem(sys, fmt="deepmd/npy")
                d.append(k)
            for k in d:
                for i in range(len(k)):
                    cell = k["cells"][i]
                    if k.nopbc:
                        cell = None
                    coord = k["coords"][i]
                    ori_atype = k["atom_types"]
                    anames = k["atom_names"]
                    atype = np.array([all_type_map.index(anames[j]) for j in ori_atype])
                    natoms = atype.shape[0]
                    e, f, v = self.evaluate(coord, cell, atype)
                    infer_energies.append(e/natoms)
                    infer_forces.extend(f.reshape(-1))
                    infer_virials.extend(v.reshape(-1)/natoms)
                    gt_energies.append(k["energies"][i]/natoms)
                    gt_forces.extend(k["forces"][i].reshape(-1))
                    gt_virials.extend(k["virials"][i].reshape(-1)/natoms)
        rmse_e = rmse(np.array(infer_energies), np.array(gt_energies)) 
        rmse_f = rmse(np.array(infer_forces), np.array(gt_forces))
        rmse_v = rmse(np.array(infer_virials), np.array(gt_virials))           
        return rmse_e, rmse_f, rmse_v

if __name__ == "__main__":
    models = glob.glob("*/model.ckpt.pt")
    valid_set = "new_data/"
    monitor_report = pd.DataFrame(columns=[
        "hparam", "steps", f"rmse_e[eV/atom]", f"rmse_f[eV/ang.]", f"rmse_v[eV/atom]", f"ava. wall_time[s/100steps]"])
    for mm in models:
        train_log = os.path.join(os.path.dirname(mm), 'train.log')
        lcurve_file = os.path.join(os.path.dirname(mm), 'lcurve.out')
        eff = get_eff(train_log)
        step_prog = get_step_prog(lcurve_file)
        if (step_prog is None) or (step_prog <1e4):
            continue
        d = DPPTPredict()
        d.load_model(mm)
        rmse_e, rmse_f, rmse_v = d.monitor_rmse(valid_set)
        hp = mm.split("/")[0]
        print(f"{hp}\t{rmse_e} {rmse_f} {rmse_v}")
        monitor_report = monitor_report._append(
            {
                "hparam": hp,
                "steps": step_prog,
                f"rmse_e[eV/atom]": rmse_e,
                f"rmse_f[eV/ang.]": rmse_f,
                f"rmse_v[eV/atom]": rmse_v,
                f"ava. wall_time[s/100steps]": eff
            },
        ignore_index=True)
    
    print(monitor_report.to_markdown())
    
