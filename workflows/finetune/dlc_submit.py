"""Generate or submit a site-specific DeepMD fine-tuning job.

The original project contained private cloud identifiers and NAS paths. This
version reads all site settings from environment variables and defaults to a
dry-run, so importing or testing it cannot submit a job accidentally.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s: %(message)s")


def parse_decaying_hparam(hparam: str) -> tuple[float | None, float | None]:
    try:
        start, end = hparam.split("->", 1)
        return float(start), float(end)
    except (ValueError, TypeError):
        return None, None


def reg_name(exp_name: str) -> str:
    return "".join(exp_name.split()).replace("->", "_to_").replace(":", "-")


def config_json(
    exp_name: str,
    template: str,
    pref_e: str,
    pref_f: str,
    pref_v: str,
    lr: str,
    n_dim: int,
    e_dim: int,
    a_dim: int,
    nlayer: int,
    bsz: str,
) -> Path:
    with open(template, encoding="utf-8") as handle:
        data = json.load(handle)
    start_e, end_e = parse_decaying_hparam(pref_e)
    start_f, end_f = parse_decaying_hparam(pref_f)
    start_v, end_v = parse_decaying_hparam(pref_v)
    start_lr, end_lr = parse_decaying_hparam(lr)
    repflow = data["model"]["descriptor"]["repflow"]
    repflow.update(n_dim=n_dim, e_dim=e_dim, a_dim=a_dim, nlayers=nlayer)
    data["loss"].update(
        start_pref_e=start_e,
        limit_pref_e=end_e,
        start_pref_f=start_f,
        limit_pref_f=end_f,
        start_pref_v=start_v,
        limit_pref_v=end_v,
    )
    data["learning_rate"].update(start_lr=start_lr, stop_lr=end_lr)
    data["training"]["training_data"]["batch_size"] = bsz
    target = Path(exp_name) / "input.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    return target


def build_command(exp_path: Path) -> list[str]:
    executable = os.environ.get("ICA_DLC_EXECUTABLE", "dlc")
    command = os.environ.get(
        "ICA_TRAIN_COMMAND",
        "dp --pt train input.json --finetune pretrain.pt > train.log 2>&1",
    )
    cmd = [executable, "submit", "pytorchjob", "--name", exp_path.name]
    options = {
        "ICA_DLC_WORKER_CPU": "--worker_cpu",
        "ICA_DLC_WORKER_GPU": "--worker_gpu",
        "ICA_DLC_WORKER_MEMORY": "--worker_memory",
        "ICA_DLC_WORKERS": "--workers",
        "ICA_DLC_IMAGE": "--worker_image",
        "ICA_DLC_DATA_SOURCES": "--data_sources",
        "ICA_DLC_WORKSPACE_ID": "--workspace_id",
        "ICA_DLC_RESOURCE_ID": "--resource_id",
    }
    for env_name, flag in options.items():
        value = os.environ.get(env_name)
        if value:
            if flag == "--worker_memory":
                value = f"{value}Gi"
            cmd.extend([flag, value])
    cmd.extend(["--command", f"cd {exp_path} && {command}"])
    return cmd


def submit(exp_path: Path, dry_run: bool = True) -> int:
    if (exp_path / "lcurve.out").exists():
        logging.info("%s already contains lcurve.out; skipping", exp_path)
        return 0
    cmd = build_command(exp_path)
    logging.info("%s", " ".join(cmd))
    if dry_run:
        return 0
    return subprocess.run(cmd, check=False).returncode


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("exp_path", type=Path)
    parser.add_argument("--submit", action="store_true", help="submit instead of dry-run")
    args = parser.parse_args()
    return submit(args.exp_path, dry_run=not args.submit)


if __name__ == "__main__":
    raise SystemExit(main())
