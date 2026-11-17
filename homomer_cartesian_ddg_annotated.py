#!/usr/bin/env python
# Copyright (c) 2026 Benjamin E. Peterson. CC BY 4.0; see LICENSE.
"""HC-ddG (Homomer Cartesian ddG): folding ddG in homo-oligomers.

This script is a PyRosetta modification of the Cartesian ddG protocol of Park
et al. (2016) and Frenz et al. (2020) for homo-oligomers. Because a
substitution in a homo-oligomer is present in every copy, the site is located
on every chain, all copies are mutated as a single set, and the assembly is
relaxed as a whole without enforcing symmetry, with each chain free to adjust
independently; the energy change is then divided by the number of chains.
Rosetta's cartesian_ddg application is not run; its routines (mutate_pose,
optimize_structure) are called from a loop written for this study. Compared
with Frenz et al., the backbone in the 8 A shell around the sites is flexible
rather than fixed, each arm runs one round instead of up to five, and the wild
type is the lowest-energy of three relaxations instead of twenty; in this
benchmark, the heavier settings did not improve accuracy and required more
time. The proline fragment sampler is disabled by default because it reduced
accuracy, rounds beyond 15 REU are counted rather than clipped, and each
substitution is reseeded; thus, results do not depend on the number of
workers. Only folding ddG is computed, not binding ddG. Usage and settings
are described in README.md and details/ADDITIONAL_DETAILS.md.
"""

import argparse
import contextlib
import csv
import json
import multiprocessing as mp
import os
import signal
import string
import sys
from pathlib import Path

import pyrosetta
from pyrosetta import Pose, create_score_function, pose_from_pdb
from pyrosetta.rosetta.core.chemical import aa_from_oneletter_code
from pyrosetta.rosetta.core.select.residue_selector import (
    NeighborhoodResidueSelector,
    ResidueIndexSelector,
)
from pyrosetta.rosetta.protocols import hybridization
from pyrosetta.rosetta.protocols.ddg.CartesianddG import (
    MutationSet,
    involves_prolines,
    mutate_pose,
    optimize_structure,
    sample_fragments,
)
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.utility import (
    vector1_core_chemical_AA,
    vector1_unsigned_long,
)

# Fixed constants: version, amino acids, the 15 REU round limit and the output
# columns.
VERSION = "1.0.0"
AA20 = "ACDEFGHIKLMNPQRSTVWY"
CLAMP = 15.0
THREAD_VARIABLES = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS")
OUT_COLUMNS = ["pdb", "mutation", "cbddg", "exp_ddg_kcal", "n_targets",
               "n_clamped"]
PER_RESIDUE_COLUMNS = ["pdb", "mutation", "resnum", "n_faces",
                       "n_faces_found", "complete", "mean_over_n_faces",
                       "face_spread", "faces"]

# Settings, read from CB_* environment variables and checked;
# details/ADDITIONAL_DETAILS.md describes each one.
STRUCT_DIR = os.environ.get("CB_STRUCT", "")
BENCH_CSV = os.environ.get("CB_BENCH", "")
OUT_CSV = os.environ.get("CB_OUT", "hcddg_out.csv")
WTDIR = Path(os.environ.get("CB_WTDIR", "hcddg_wt_cache"))
SFXN = os.environ.get("CB_SFXN", "ref2015_cart")
PACK_RADIUS = float(os.environ.get("CB_PACK", "6.0"))
BBNBRS = int(os.environ.get("CB_BBNBRS", "1"))
ROUNDS = int(os.environ.get("CB_ROUNDS", "1"))
WT_REPEATS = int(os.environ.get("CB_WTREPEATS", "5"))
NPREP = int(os.environ.get("CB_NPREP", "3"))
PREP_CART = os.environ.get("CB_PREPCART", "1") == "1"
PREP_RESTRAIN = os.environ.get("CB_PREP_RESTRAIN", "none").lower()
RESEED_BASE = int(os.environ.get("CB_RESEED_BASE", "1000"))
PREP_SEED = int(os.environ.get("CB_PREP_SEED", "7"))
NFRAGS = int(os.environ.get("CB_NFRAGS", "30"))
FRAG_NBRS = int(os.environ.get("CB_FRAG_NBRS", "2"))
NTRIALS = int(os.environ.get("CB_NTRIALS", "30"))
EXPECT_N = int(os.environ.get("CB_EXPECT_N", "0"))
SHELL = os.environ.get("CB_SHELL", "selector8").lower()
SHELL_RADIUS = float(os.environ.get("CB_SHELL_RADIUS", "8.0"))

BOOLEAN_WORDS = {"1": True, "true": True, "yes": True, "on": True,
                 "0": False, "false": False, "no": False, "off": False}
BOOLEAN_DEFAULTS = {"CB_RESEED": "1", "CB_OPTIMIZE_PROLINE": "0",
                    "CB_PER_RESIDUE": "0", "CB_MISSING_DENSITY": "0",
                    "CB_ZERO_OCC": "0", "CB_IGNORE_WATERS": "1",
                    "CB_FLEXBB": "1", "CB_CLAMP": "0"}
BOOLEANS = {name: BOOLEAN_WORDS.get(
                os.environ.get(name, default).strip().lower())
            for name, default in BOOLEAN_DEFAULTS.items()}
BAD_BOOLEANS = [name for name, value in BOOLEANS.items() if value is None]
if BAD_BOOLEANS:
    raise SystemExit("Not a boolean: "
                     + ", ".join(f"{name}={os.environ[name]!r}"
                                 for name in BAD_BOOLEANS)
                     + ". Use 0/1, true/false, yes/no or on/off.")
RESEED = BOOLEANS["CB_RESEED"]
OPTIMIZE_PROLINE = BOOLEANS["CB_OPTIMIZE_PROLINE"]
PER_RESIDUE = BOOLEANS["CB_PER_RESIDUE"]
MISSING_DENSITY = BOOLEANS["CB_MISSING_DENSITY"]
ZERO_OCC = BOOLEANS["CB_ZERO_OCC"]
IGNORE_WATERS = BOOLEANS["CB_IGNORE_WATERS"]
FLEXBB = BOOLEANS["CB_FLEXBB"]
CLAMP_APPLY = BOOLEANS["CB_CLAMP"]

if not 1 <= NPREP <= 20:
    raise SystemExit(f"CB_NPREP={NPREP} refused; allowed range is 1 to 20 "
                     "(Frenz uses 20).")
if not 0 <= BBNBRS <= 5:
    raise SystemExit(f"CB_BBNBRS={BBNBRS} refused; allowed range is 0 to 5.")
if PREP_RESTRAIN not in ("none", "bb", "bbsc"):
    raise SystemExit(f"CB_PREP_RESTRAIN={PREP_RESTRAIN!r}; use none, bb or "
                     "bbsc.")
if SHELL not in ("nbrradius6", "selector8"):
    raise SystemExit(f"CB_SHELL={SHELL!r}; use nbrradius6 or selector8.")
if SHELL == "selector8" and not 3.0 <= SHELL_RADIUS <= 20.0:
    raise SystemExit(f"CB_SHELL_RADIUS={SHELL_RADIUS} refused; use 3 to "
                     "20 A.")

# One-line summary of the settings, printed with each run.
if SHELL == "selector8":
    SHELL_TAG = "@" + str(SHELL_RADIUS)
else:
    SHELL_TAG = "+pack" + str(PACK_RADIUS)
PROTOCOL_TAG = (
    f"v{VERSION} prep[n{NPREP},r{WT_REPEATS},"
    f"{'cart' if PREP_CART else 'tors'},restrain={PREP_RESTRAIN}] "
    f"ddg[rounds={ROUNDS},flex_bb={FLEXBB},shell={SHELL}{SHELL_TAG},"
    f"bbnbrs={BBNBRS},clamp={'on' if CLAMP_APPLY else 'count-only'},"
    f"proline={'on' if OPTIMIZE_PROLINE else 'off'}] "
    f"reseed_base={RESEED_BASE} "
    f"reseed={'per-mutation' if RESEED else 'per-worker'}"
)


# Rosetta start-up options shared by every worker process.
def _init_opts():
    beta = " -corrections::beta_nov16 true" if "beta_nov16" in SFXN else ""
    return ("-ex1 -ex2 -use_input_sc -flip_HNQ -no_optH false "
            "-ignore_unrecognized_res "
            f"-mute all -fa_max_dis 9.0 -score:weights {SFXN} "
            f"{'-missing_density_to_jump ' if MISSING_DENSITY else ''}"
            f"{'-ignore_zero_occupancy false ' if ZERO_OCC else ''}"
            f"{'-ignore_waters ' if IGNORE_WATERS else ''}"
            "-relax:default_repeats 5 "
            f"-constant_seed -jran 7{beta}")


def _init():
    pyrosetta.init(extra_options=_init_opts())


# Wild-type preparation: relax each structure NPREP times, keep the
# lowest-energy model and cache it.
def _wt_cache_path(pdb):
    seed_tag = f"p{PREP_SEED}" if RESEED else "pfree"
    import_tag = (f"md{int(MISSING_DENSITY)}zo{int(ZERO_OCC)}"
                  f"iw{int(IGNORE_WATERS)}")
    tag = (f"{SFXN}_r{WT_REPEATS}_n{NPREP}_"
           f"{'cart' if PREP_CART else 'tors'}"
           f"_restrain-{PREP_RESTRAIN}_{seed_tag}_{import_tag}")
    return str(WTDIR / f"{pdb}_wt__{tag}.pdb")


def prep_wt(pdb):
    WTDIR.mkdir(parents=True, exist_ok=True)
    wt_path = _wt_cache_path(pdb)
    source_path = wt_path + ".source"
    if os.path.exists(wt_path) and os.path.exists(source_path):
        with open(source_path) as source_file:
            previous = source_file.read().strip()
        if previous != str(STRUCT_DIR):
            raise RuntimeError(
                f"WT cache collision for {pdb}: cached pose was prepared "
                f"from\n  {previous}\nbut CB_STRUCT is now\n  {STRUCT_DIR}\n"
                "Same filename, different input structure. Use a different "
                "CB_WTDIR.")
    if os.path.exists(wt_path):
        return (pdb, wt_path)
    if RESEED:
        pyrosetta.rosetta.basic.random.init_random_generators(PREP_SEED,
                                                              "mt19937")
    score_function = create_score_function(SFXN)
    pose = pose_from_pdb(f"{STRUCT_DIR}/{pdb}.pdb")
    best_pose, best_score = None, float("inf")
    for _ in range(NPREP):
        trial = Pose()
        trial.assign(pose)
        relax_function = score_function.clone()
        relax = FastRelax(relax_function, WT_REPEATS)
        relax.cartesian(PREP_CART)
        if PREP_RESTRAIN in ("bb", "bbsc"):
            relax.constrain_coords(True)
            relax.ramp_down_constraints(False)
            if PREP_RESTRAIN == "bbsc":
                relax.coord_constrain_sidechains(True)
        relax.apply(trial)
        score = float(score_function(trial))
        if score < best_score:
            best_pose, best_score = trial, score
    best_pose.dump_pdb(wt_path)
    with open(source_path, "w") as source_file:
        source_file.write(str(STRUCT_DIR))
    return (pdb, wt_path)


# Locate the site on every chain and build the shell of residues optimised
# around it.
def parse_resnum(resnum):
    text = str(resnum).strip()
    number = text.rstrip(string.ascii_letters)
    icode = text[len(number):]
    if not number.lstrip("-").isdigit() or len(icode) > 1:
        raise ValueError(f"resnum {resnum!r} is not a residue number with an "
                         "optional one-letter insertion code")
    return int(number), (icode or " ")


def resnum_order(resnum):
    number, icode = parse_resnum(resnum)
    return number, icode.strip()


def both_positions(pose, resnum, wt):
    number, icode = parse_resnum(resnum)
    info = pose.pdb_info()
    chains = {info.chain(index)
              for index in range(1, pose.total_residue() + 1)}
    indices = [info.pdb2pose(chain, number, icode) for chain in chains]
    return sorted(index for index in indices
                  if index != 0 and pose.residue(index).name1() == wt)


def build_shell(pose, targets, threshold):
    shell = vector1_unsigned_long()
    if SHELL == "selector8":
        focus = ResidueIndexSelector(",".join(str(target)
                                              for target in targets))
        selector = NeighborhoodResidueSelector(focus, SHELL_RADIUS, True)
        mask = selector.apply(pose)
        for index in range(1, pose.total_residue() + 1):
            if mask[index]:
                shell.append(index)
        return shell
    target_set = set(targets)
    for index in range(1, pose.total_residue() + 1):
        if index in target_set:
            continue
        residue = pose.residue(index)
        for target in targets:
            target_residue = pose.residue(target)
            cutoff = (residue.nbr_radius() + target_residue.nbr_radius()
                      + threshold)
            distance_squared = residue.nbr_atom_xyz().distance_squared(
                target_residue.nbr_atom_xyz())
            if distance_squared <= cutoff * cutoff:
                shell.append(index)
                break
    return shell


# Optional backbone fragments for proline substitutions (off by default).
def attach_proline_fragments(pose, mutation_set, targets, mut):
    prolines = involves_prolines(pose, mutation_set)
    if not prolines:
        return 0
    letters = list(pose.sequence())
    for target in targets:
        letters[target - 1] = mut
    sequence = "".join(letters)
    n_decorated = 0
    n_residues = pose.total_residue()
    for position in prolines:
        lower = max(1, int(position) - FRAG_NBRS)
        upper = min(n_residues, int(position) + FRAG_NBRS)
        fragment_size = upper - lower + 1
        if fragment_size < 4:
            continue
        fragment_set = hybridization.create_fragment_set_no_ssbias(
            sequence[lower - 1:upper], fragment_size, NFRAGS, "D")
        mutation_set.add_fragments(lower, fragment_set)
        n_decorated += 1
    return n_decorated


# Score substitutions: for each one, optimise a wild-type arm and a mutant arm
# over the shell and report the energy difference per protomer.
def score_chunk(task):
    pdb, wt_path, substitutions, mode = task
    score_function = create_score_function(SFXN)
    wt_pose = pose_from_pdb(wt_path)
    rows = []
    per_residue_rows = []
    for label, wt, resnum, mut, exp_value in substitutions:
        if RESEED:
            pyrosetta.rosetta.numeric.random.rg().set_seed("mt19937",
                                                           RESEED_BASE)
        # Find the site on every chain (the first chain only in single mode),
        # check the count and build the shell.
        found = both_positions(wt_pose, resnum, wt)
        if not found:
            print(f"  [skip] {pdb} {label}: no chain has {wt} at residue "
                  f"{resnum}.", flush=True)
            continue
        if EXPECT_N != 0 and len(found) != EXPECT_N:
            print(f"  [FAIL] {pdb} {label}: found {len(found)} protomer "
                  f"target(s), expected {EXPECT_N}. Row dropped, not "
                  "averaged.", flush=True)
            continue
        if EXPECT_N == 0 and len(found) < 2:
            print(f"  [WARN] {pdb} {label}: found in {len(found)} protomer "
                  "only. Set CB_EXPECT_N to require every protomer.",
                  flush=True)
        targets = found[:1] if mode == "single" else found
        div = float(len(targets)) if mode in ("homomer", "bilateral") else 1.0
        shell = build_shell(wt_pose, targets, PACK_RADIUS)
        target_indices = vector1_unsigned_long()
        for index in targets:
            target_indices.append(index)
        n_clamped = 0
        round_values = []
        for _ in range(ROUNDS):
            # Optimise the wild-type arm, then mutate every target and optimise
            # the mutant arm.
            wt_types = vector1_core_chemical_AA()
            for _ in targets:
                wt_types.append(aa_from_oneletter_code(wt))
            wt_set = MutationSet(target_indices, wt_types, 1)
            wt_arm = Pose()
            wt_arm.assign(wt_pose)
            optimize_structure(wt_set, score_function, wt_arm, shell, FLEXBB,
                               True, BBNBRS)
            wt_energy = float(score_function(wt_arm))

            mut_types = vector1_core_chemical_AA()
            for _ in targets:
                mut_types.append(aa_from_oneletter_code(mut))
            mut_set = MutationSet(target_indices, mut_types, 1)
            n_prolines = 0
            if OPTIMIZE_PROLINE:
                n_prolines = attach_proline_fragments(wt_pose, mut_set,
                                                      targets, mut)
            mut_arm = Pose()
            mut_arm.assign(wt_pose)
            mutate_pose(mut_arm, mut_set, score_function)
            if n_prolines > 0:
                sample_fragments(mut_arm, mut_set, score_function, BBNBRS,
                                 NTRIALS)
            optimize_structure(mut_set, score_function, mut_arm, shell,
                               FLEXBB, True, BBNBRS)
            mut_energy = float(score_function(mut_arm))

            # ddG per protomer for this round; values beyond 15 REU are
            # counted, and clipped only with CB_CLAMP=1.
            raw = (mut_energy - wt_energy) / div
            if abs(raw) > CLAMP:
                n_clamped += 1
            round_values.append(max(-CLAMP, min(CLAMP, raw))
                                if CLAMP_APPLY else raw)

            # Optional energy change of each shell residue on each protomer.
            if PER_RESIDUE:
                info = wt_pose.pdb_info()
                by_residue = {}
                for index in shell:
                    delta = float(
                        mut_arm.energies().residue_total_energy(index)
                        - wt_arm.energies().residue_total_energy(index))
                    key = (int(info.number(index)),
                           info.icode(index).strip())
                    by_residue.setdefault(key, {})[info.chain(index)] = (
                        round(delta, 4))
                n_faces = len(targets)
                for (number, icode), faces in sorted(by_residue.items()):
                    complete = len(faces) == n_faces
                    mean = round(sum(faces.values()) / len(faces), 4)
                    spread = round(max(faces.values()) - min(faces.values()),
                                   4)
                    residue_label = (f"{number}{icode}" if icode != ""
                                     else number)
                    per_residue_rows.append((
                        pdb, label, residue_label, n_faces, len(faces),
                        "yes" if complete else "no",
                        mean if complete else "", spread,
                        ";".join(f"{chain}={value}"
                                 for chain, value in sorted(faces.items()))))
        mean_value = round(sum(round_values) / len(round_values), 4)
        rows.append((pdb, label, mean_value, exp_value, len(targets),
                     n_clamped))
    return rows, per_residue_rows


# Report and skip a chunk that fails instead of stopping the run.
def safe_score_chunk(task):
    try:
        return score_chunk(task)
    except Exception as error:
        print(f"  [WARN] chunk {task[0]} ({len(task[2])} muts) FAILED, "
              f"skipped: {type(error).__name__}: {str(error)[:100]}",
              flush=True)
        return [], []


# Record the settings in force in <out>.protocol.json.
def write_sidecar(out, mode):
    settings = {
        "version": VERSION,
        "tag": PROTOCOL_TAG,
        "sfxn": SFXN,
        "nprep": NPREP,
        "wt_repeats": WT_REPEATS,
        "rounds": ROUNDS,
        "flex_bb": FLEXBB,
        "shell": SHELL,
        "shell_radius": SHELL_RADIUS if SHELL == "selector8" else None,
        "pack_radius": PACK_RADIUS if SHELL != "selector8" else None,
        "struct_dir": str(STRUCT_DIR),
        "bbnbrs": BBNBRS,
        "prep_restrain": PREP_RESTRAIN,
        "optimize_proline": OPTIMIZE_PROLINE,
        "zero_occ": ZERO_OCC,
        "ignore_waters": IGNORE_WATERS,
        "missing_density": MISSING_DENSITY,
        "expect_n": EXPECT_N,
        "prep_cartesian": PREP_CART,
        "reseed": RESEED,
        "reseed_base": RESEED_BASE,
        "clamp_reu": CLAMP,
        "prep_seed": PREP_SEED,
        "clamp_applied": CLAMP_APPLY,
        "mode": mode,
    }
    try:
        with open(str(out) + ".protocol.json", "w") as sidecar_file:
            json.dump(settings, sidecar_file, indent=2)
    except OSError as error:
        print(f"  [WARN] could not write protocol sidecar: {error}",
              flush=True)


# Command line: read the substitutions, prepare the wild types, score the
# substitutions in parallel and write the output.
def main():
    global STRUCT_DIR
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("structures", nargs="?", default=None,
                        help="directory of input structures "
                             "(default: CB_STRUCT)")
    parser.add_argument("mutations", nargs="?", default=None,
                        help="substitution CSV (default: --bench or "
                             "CB_BENCH)")
    parser.add_argument("--mode", default="homomer",
                        choices=["homomer", "bilateral", "single"],
                        help="homomer: mutate every protomer, divide by n "
                             "(bilateral is an alias). single: mutate one "
                             "protomer only, the control.")
    parser.add_argument("--pdbs", default=None,
                        help="comma-separated PDB IDs to score "
                             "(default: all)")
    parser.add_argument("--bench", default=None,
                        help="mutation CSV (default: CB_BENCH)")
    parser.add_argument("--out", default=None,
                        help="output CSV (default: CB_OUT or hcddg_out.csv)")
    parser.add_argument("--workers", type=int,
                        default=max(1, (os.cpu_count() or 1) - 2),
                        help="parallel worker processes (default: CPUs - 2)")
    args = parser.parse_args()
    signal.signal(signal.SIGTERM,
                  lambda signum, frame: sys.exit(128 + signum))
    for variable in THREAD_VARIABLES:
        os.environ[variable] = "1"

    # Take the inputs from the command line or the CB_* variables, read the
    # substitution CSV and group the rows by structure.
    if args.structures is not None:
        STRUCT_DIR = args.structures
        os.environ["CB_STRUCT"] = STRUCT_DIR
    bench = args.mutations or args.bench or BENCH_CSV
    out = args.out or OUT_CSV
    if not STRUCT_DIR or not bench:
        raise SystemExit("Give the structure directory and the substitution "
                         "CSV, for example: python homomer_cartesian_ddg.py "
                         "structures mutations.csv")
    selected = args.pdbs.split(",") if args.pdbs is not None else []
    with open(bench, newline="") as bench_file:
        reader = csv.DictReader(bench_file)
        missing = [column
                   for column in ("pdb", "mutation", "wt", "resnum", "mut")
                   if column not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(f"{bench} lacks column(s): {', '.join(missing)}")
        bench_rows = [row for row in reader
                      if not selected or row["pdb"] in selected]
    by_pdb = {}
    for row in bench_rows:
        by_pdb.setdefault(row["pdb"], []).append(row)
    if not by_pdb:
        raise SystemExit(
            f"No rows in {bench} matched "
            + (f"--pdbs {args.pdbs}" if selected else "the benchmark")
            + ". Check the PDB IDs against the pdb column.")

    # Prepare each wild type once, in parallel.
    mp.set_start_method("spawn", force=True)
    with mp.Pool(min(len(by_pdb), args.workers), initializer=_init) as pool:
        wt_paths = dict(pool.map(prep_wt, list(by_pdb)))
    print("WT prep done:", list(wt_paths), flush=True)

    # Check each substitution and split them into chunks of four for the
    # workers.
    tasks = []
    for pdb, pdb_rows in by_pdb.items():
        if pdb not in wt_paths:
            continue
        substitutions = []
        for row in pdb_rows:
            wt = row["wt"].strip().upper()
            mut = row["mut"].strip().upper()
            try:
                parse_resnum(row["resnum"])
            except ValueError as error:
                print(f"  [skip] {pdb} {row['mutation']}: {error}",
                      flush=True)
                continue
            if (len(wt) != 1 or len(mut) != 1 or wt not in AA20
                    or mut not in AA20):
                print(f"  [skip] {pdb} {row['mutation']}: wt {row['wt']!r} "
                      f"or mut {row['mut']!r} is not a standard one-letter "
                      "amino acid", flush=True)
                continue
            exp_text = (row.get("exp_ddg_kcal") or "").strip()
            exp_value = ""
            if exp_text != "":
                try:
                    exp_value = float(exp_text)
                except ValueError:
                    print(f"  [note] {pdb} {row['mutation']}: exp_ddg_kcal "
                          f"{exp_text!r} is not a number; written blank.",
                          flush=True)
            substitutions.append((row["mutation"], wt, row["resnum"], mut,
                                  exp_value))
        for start in range(0, len(substitutions), 4):
            tasks.append((pdb, wt_paths[pdb], substitutions[start:start + 4],
                          args.mode))
    print(f"{len(tasks)} chunks; homomer cartesian ddG, {ROUNDS} rounds/arm, "
          f"WT prep {NPREP}x{WT_REPEATS} repeats", flush=True)
    print(f"PROTOCOL: {PROTOCOL_TAG}", flush=True)

    # Score the chunks in parallel and write the output sorted, even if the run
    # is interrupted.
    done = 0
    written = 0
    results = []
    per_residue_results = []
    with contextlib.ExitStack() as stack:
        out_file = stack.enter_context(open(out, "w", newline=""))
        writer = csv.writer(out_file)
        writer.writerow(OUT_COLUMNS)
        per_residue_file = None
        per_residue_writer = None
        if PER_RESIDUE:
            per_residue_file = stack.enter_context(
                open(str(out) + ".perresidue.csv", "w", newline=""))
            per_residue_writer = csv.writer(per_residue_file)
            per_residue_writer.writerow(PER_RESIDUE_COLUMNS)
        try:
            with mp.Pool(args.workers, initializer=_init,
                         maxtasksperchild=1) as pool:
                for chunk_rows, chunk_residue_rows in pool.imap_unordered(
                        safe_score_chunk, tasks):
                    results.extend(chunk_rows)
                    per_residue_results.extend(chunk_residue_rows)
                    written += len(chunk_rows)
                    done += 1
                    print(f"  [{done}/{len(tasks)} chunks]  "
                          f"+{len(chunk_rows)} rows", flush=True)
        finally:
            for row in sorted(results, key=lambda row: (
                    str(row[0]), str(row[1]), tuple(map(str, row)))):
                writer.writerow(row)
            out_file.flush()
            if per_residue_writer is not None:
                for row in sorted(per_residue_results, key=lambda row: (
                        str(row[0]), str(row[1]), resnum_order(row[2]),
                        tuple(map(str, row)))):
                    per_residue_writer.writerow(row)
                per_residue_file.flush()
            write_sidecar(out, args.mode)

    # Warn if any substitution produced no result.
    expected = sum(len(task[2]) for task in tasks)
    if written != expected:
        print(f"\n[warning] {expected - written} of {expected} mutations "
              "produced no result. Look for [WARN], [FAIL] or [skip] lines "
              "above before using this table.", flush=True)


if __name__ == "__main__":
    main()
