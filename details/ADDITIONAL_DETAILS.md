# HC-ddG: additional details

## Input

**Structures.** One PDB file per protein, named `<ID>.pdb` (for example `1ROP.pdb`), in one folder. Each copy of the chain must be a separate chain in the file, with the same residue numbering in every chain. Remove any chains that are not copies.

**Substitutions.** A CSV file with one row per substitution, for example `example/mutations.csv`:

    pdb,mutation,wt,resnum,mut
    1ROP,L41V,L,41,V
    1ROP,L41A,L,41,A
    1ROP,D30G,D,30,G

| Column | Meaning |
|---|---|
| `pdb` | Structure file name without `.pdb` |
| `mutation` | Any label; it is copied to the output |
| `wt` | Original amino acid (one letter) |
| `resnum` | Residue number in the PDB file, with an optional insertion code such as `52A` |
| `mut` | New amino acid (one letter) |
| `exp_ddg_kcal` | Optional measured ddG; it is copied to the output |

## Running

python homomer_cartesian_ddg.py structures mutations.csv

`structures` is the folder of structure files and `mutations.csv` is the substitution CSV. The results are written to `hcddg_out.csv`.

## Output

`hcddg_out.csv` has one row per substitution:

| Column | Meaning |
|---|---|
| `pdb`, `mutation` | Copied from the input |
| `cbddg` | Predicted ddG per chain, in Rosetta energy units (REU). Positive values mean the substitution is destabilizing |
| `exp_ddg_kcal` | Measured ddG, if it was given |
| `n_targets` | Number of chains that were mutated |
| `n_clamped` | Number of rounds with a value beyond ±15 REU |

The results are in Rosetta energy units, not kcal/mol. They are meant for ranking substitutions within one protein.

Two more files are written next to it: `hcddg_out.csv.protocol.json`, which records the settings used, and, if `CB_PER_RESIDUE=1` is set, `hcddg_out.csv.perresidue.csv`, which gives the energy change of each nearby residue on each chain. The output is saved even if the run is stopped early (Ctrl-C).

## Options

| Option | What it does |
|---|---|
| `--out FILE` | Output file name (default `hcddg_out.csv`) |
| `--pdbs 1ROP,2ZTA` | Score only these structures from the CSV (default: all) |
| `--workers N` | Number of parallel processes (default: CPUs minus two) |
| `--mode single` | Mutate only one chain, as a control (default: every chain) |

To skip any substitution whose site is not found on every chain, give the number of chains before the command:

    CB_EXPECT_N=2 python homomer_cartesian_ddg.py structures mutations.csv

Older ways of giving the inputs (`CB_STRUCT`, `CB_BENCH` and `--bench`) still work.

## Advanced settings

Most users do not need these. They are set as environment variables, in the same way as `CB_EXPECT_N` above. Defaults are in brackets.

| Variable | Meaning |
|---|---|
| `CB_OUT` [`hcddg_out.csv`] | Output file name, if `--out` is not given |
| `CB_EXPECT_N` [0] | Number of chains that must carry each site; rows that do not match are skipped. 0 turns the check off |
| `CB_WTDIR` [`hcddg_wt_cache`] | Folder where prepared starting structures are cached |
| `CB_PREP_SEED` [7] | Random seed for preparing the starting structure |
| `CB_NPREP` [3] | Number of relaxations of the starting structure; the lowest-energy one is kept |
| `CB_WTREPEATS` [5] | FastRelax repeats in each of those relaxations |
| `CB_PREPCART` [1] | 1 relaxes the starting structure in Cartesian space; any other value uses torsion space |
| `CB_PREP_RESTRAIN` [`none`] | Coordinate restraints while preparing the starting structure: `none`, `bb` or `bbsc` |
| `CB_RESEED_BASE` [1000] | Random seed set before every substitution |
| `CB_RESEED` [1] | 0 turns off both random seeds |
| `CB_ROUNDS` [1] | Rounds per model, averaged |
| `CB_SHELL` [`selector8`] | How nearby residues are chosen: `selector8` (all residues within a radius) or `nbrradius6` (Rosetta's own rule) |
| `CB_SHELL_RADIUS` [8.0] | Radius for `selector8`, in Å (3 to 20) |
| `CB_PACK` [6.0] | Distance used by `nbrradius6` |
| `CB_FLEXBB` [1] | 1 lets the backbone of nearby residues move |
| `CB_BBNBRS` [1] | Backbone residues freed on each side of a site (0 to 5) |
| `CB_SFXN` [`ref2015_cart`] | Rosetta score function |
| `CB_CLAMP` [0] | 1 caps each round at ±15 REU; otherwise those rounds are only counted |
| `CB_OPTIMIZE_PROLINE` [0] | 1 turns on extra backbone sampling for proline substitutions (with `CB_NFRAGS` [30], `CB_FRAG_NBRS` [2] and `CB_NTRIALS` [30]) |
| `CB_PER_RESIDUE` [0] | 1 also writes the per-residue, per-chain energy file |
| `CB_MISSING_DENSITY` [0] | Sets Rosetta's `-missing_density_to_jump` |
| `CB_ZERO_OCC` [0] | 1 keeps atoms with zero occupancy |
| `CB_IGNORE_WATERS` [1] | 1 ignores water molecules |

## Differences from the standard protocols

HC-ddG follows the Cartesian ddG protocol of Park et al. (2016) and Frenz et al. (2020), as run by Rosetta's `cartesian_ddg` application. The older protocol of Kellogg et al. (2011), run by `ddg_monomer`, differs more: it restrains C-alpha distances and runs 50 iterations per arm. Both applications are documented for monomers. HC-ddG changes the following.

- The substitution is made in every chain that carries the wild-type residue at the site, as one mutation set, and the energy difference is divided by n. No symmetry is imposed. `cartesian_ddg` can take several sites from a mutation file, but the user lists them. HC-ddG finds the site on each chain and, when `CB_EXPECT_N` is set, checks the count.
- The outer loop is our own. `optimize_structure` and `mutate_pose` are called directly, and Rosetta's `cartesian_ddg` application is not run.
- Folding ddG only. The interface mode of `cartesian_ddg`, which its documentation says is not thoroughly benchmarked, is not included.
- Input is a CSV of substitutions rather than a mutation file. Output is one CSV row per substitution, with a `.protocol.json` record of the settings. With `CB_PER_RESIDUE=1` the script also writes the energy change of each shell residue on each protomer, which is not part of the `cartesian_ddg` output.
- Compared with Kellogg et al., there are no C-alpha distance restraints, and each arm runs one round instead of 50 iterations scored by the mean of the top three models.
- Compared with Frenz et al.:
  - The shell is a `NeighborhoodResidueSelector` at 8 Å that includes the substituted sites. Rosetta's rule adds the neighbour radii of both residues and 6 Å, and leaves the sites out of the neighbour list.
  - The backbone is flexible in the shell, each arm runs one round, and the wild type is the lowest of three relaxations. Frenz et al. hold the backbone fixed, run up to five rounds (stopping early when the two lowest energies agree within 1 REU) and take the lowest of 20 wild-type relaxations. In our benchmark those choices did not improve accuracy.
  - `-fa_max_dis 9.0` is set, as in the Frenz example in the Rosetta documentation.
  - Rounds beyond ±15 REU are counted in `n_clamped` and left unclipped, unless `CB_CLAMP=1`.
  - The proline fragment sampler is off by default because it lowered accuracy.
  - Each substitution is reseeded and the wild-type preparation uses a fixed seed, so results do not depend on the number of workers.

The nearest equivalent of the Frenz et al. settings is `CB_ROUNDS=5 CB_WTREPEATS=4 CB_NPREP=20 CB_FLEXBB=0 CB_SHELL=nbrradius6 CB_MISSING_DENSITY=1 CB_ZERO_OCC=1 CB_IGNORE_WATERS=0`.

## Limitations

- Every chain that has the wild-type residue at the given number is mutated. Supply only copies of the protomer and set `CB_EXPECT_N`, so that a site found on the wrong number of chains is dropped.
- The protomers are relaxed together, with one shell, one MoveMap and one FastRelax run per arm. Per-protomer values in the per-residue output describe that single joint relaxation and are not independent repeats.
- The prepared wild type has not had any symmetry enforced, so relaxing without restraints lets identical chains differ. Which protomer scores higher changes with `CB_PREP_SEED`, so only the size of a difference carries information.
- For many substitutions the ddG stage returns the same value whatever its seed. Changing `CB_PREP_SEED` moves results, so use several preparation seeds for an ensemble.
- Residues Rosetta cannot read are dropped and waters are ignored by default. Ligands and metal ions that are read become part of the shell. Check structures with HETATM records.
- Values are in REU and are compressed relative to experiment. Use them to rank substitutions, or recalibrate before converting to kcal/mol.
- Values depend on the PyRosetta build. Compare results only within one build. The expected output in `example/` was produced with 2026.21+release.36f0c659eb.
- The wild-type cache is keyed on the structure file name and directory, not on the file contents. Use a new `CB_WTDIR` after editing a structure.
- HC-ddG computes folding ddG for assemblies of two to four protomers. It does not compute binding ddG.
- The systematic review mentioned in the abstract is in progress. Its screening has not been audited by the authors and its results are preliminary. The program and the benchmark do not depend on it.
