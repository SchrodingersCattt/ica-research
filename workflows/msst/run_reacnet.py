#!/usr/bin/env python3
"""Run the public ReacNetGenerator command for one LAMMPS trajectory.

The wrapper keeps all paths supplied by the caller and writes a path-free
``run_manifest.json``.  It deliberately does not ship a trajectory or a
ReacNetGenerator executable.  The command defaults mirror the manuscript:
ASE-based bond perception, Cu--N connectivity excluded, and no HMM filter.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


PERIODIC_TABLE = (
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu "
    "Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs "
    "Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl "
    "Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh "
    "Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_command(args: argparse.Namespace, input_name: str) -> list[str]:
    command = [
        args.executable,
        "-i",
        input_name,
        "--type",
        "dump",
        "-a",
        *PERIODIC_TABLE,
        "-n",
        str(args.nproc),
        "--nohmm",
        "--use-ase",
        "--ase-cutoff-mult",
        str(args.ase_cutoff_mult),
        "--ase-pair-cutoffs",
        args.ase_pair_cutoffs,
        "--stepinterval",
        str(args.stepinterval),
        "--maxspecies",
        str(args.maxspecies),
        "--items",
        args.items,
    ]
    return command


def make_link(source: Path, directory: Path) -> Path:
    """Materialize a same-volume link without copying a large trajectory."""

    target = directory / source.name
    try:
        os.link(source, target)
    except OSError:
        try:
            target.symlink_to(source)
        except OSError:
            shutil.copy2(source, target)
    return target


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--executable", default="reacnetgenerator")
    parser.add_argument("--nproc", type=int, default=1)
    parser.add_argument("--stepinterval", type=int, default=1)
    parser.add_argument("--maxspecies", type=int, default=50)
    parser.add_argument("--ase-cutoff-mult", type=float, default=1.2)
    parser.add_argument(
        "--ase-pair-cutoffs",
        default="Cu-N:0",
        help="ASE pair-cutoff overrides; Cu-N:0 excludes Cu--N connectivity.",
    )
    parser.add_argument(
        "--items",
        default="species,reactions,network,report",
        help="Comma-separated ReacNetGenerator output items.",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.nproc < 1 or args.stepinterval < 1 or args.maxspecies < 1:
        parser.error("nproc, stepinterval and maxspecies must be positive")
    return args


def main() -> int:
    args = parse_args()
    trajectory = args.trajectory.expanduser().resolve()
    if not trajectory.is_file() or trajectory.stat().st_size == 0:
        raise FileNotFoundError(f"trajectory is missing or empty: {trajectory}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    command = build_command(args, trajectory.name)
    if args.dry_run:
        print(json.dumps(command, indent=2))
        return 0

    with tempfile.TemporaryDirectory(prefix=".reacnet-", dir=args.output_dir) as work:
        workdir = Path(work)
        linked = make_link(trajectory, workdir)
        completed = subprocess.run(
            command,
            cwd=workdir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        (args.output_dir / "reacnet.log").write_text(completed.stdout, encoding="utf-8")
        if completed.returncode:
            raise RuntimeError(
                f"ReacNetGenerator failed with return code {completed.returncode}; "
                "see reacnet.log"
            )

        generated = []
        for candidate in workdir.iterdir():
            if candidate == linked or not candidate.is_file():
                continue
            target = args.output_dir / candidate.name
            shutil.move(str(candidate), str(target))
            generated.append(target.name)

    manifest_command = list(command)
    manifest_command[0] = Path(manifest_command[0]).name
    manifest = {
        "trajectory": trajectory.name,
        "trajectory_sha256": sha256(trajectory),
        "command": manifest_command,
        "parameters": {
            "nproc": args.nproc,
            "stepinterval": args.stepinterval,
            "maxspecies": args.maxspecies,
            "ase_cutoff_mult": args.ase_cutoff_mult,
            "ase_pair_cutoffs": args.ase_pair_cutoffs,
            "items": args.items,
            "hmm": False,
        },
        "generated_files": sorted(generated),
    }
    (args.output_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"output_dir": args.output_dir.name, "files": sorted(generated)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
