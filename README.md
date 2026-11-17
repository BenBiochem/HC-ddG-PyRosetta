# HC-ddG

HC-ddG (Homomer Cartesian ddG) estimates the change in folding stability caused by a point substitution in a homo-oligomer. It is a Python script that runs Rosetta's Cartesian ddG routine through PyRosetta. The substitution is made in every chain and the protomers are relaxed without imposing symmetry.

## Abstract

Many proteins function as homo-oligomers, assemblies of identical chains encoded by a single gene; thus, a substitution is present in every copy. Identical sequences need not give identical structures. Native protomers can adopt slightly different functional states, so the predicted change in folding stability (ddG) of a substitution may depend on which protomer is mutated and how the copies are relaxed. A preliminary systematic review of the literature reveals repeated evidence of residue-level differences between protomers in structure, dynamics, and function, although global asymmetry remains uncommon. Most computational ddG predictors are developed and benchmarked on monomers, where there are no copies to differ or interact. Structure-based workflows treat a homo-oligomer in one of three ways: mutating one chain, which models a mixed assembly of wild-type and mutant chains; mutating every chain; or enforcing symmetry among the copies. To our knowledge, these approaches have not been systematically compared on the same set of substitutions. This study compares these strategies on 20 homo-oligomers with 420 measured substitutions and introduces HC-ddG (Homomer Cartesian ddG), a PyRosetta modification of the Frenz et al. protocol that incorporates these findings. In Rosetta, mutating every chain and relaxing the protomers independently within a single assembly provided the most accurate ranking of substitutions, with improvements observed both near and distant from the interface. Enforcing symmetry, where feasible, resulted in the greatest loss of accuracy. HC-ddG performed comparably to a faithful reimplementation of the original Frenz et al. protocol, which also mutated every chain, and outperformed ThermoMPNN. These results show that homo-oligomer ddG is most accurately predicted by mutating every chain and relaxing the assembly as a whole, with each chain free to adjust independently. Together with the review, this indicates that protomers can occupy different energetic states at the residue level and that these differences must be considered when calculating ddG.

## How it works

The substitution is made in every protomer at once, as one mutation set, and the assembly is relaxed without imposing symmetry. The reported value is per protomer:

    ddG = (E_mut - E_wt) / n

in Rosetta energy units (REU), where n is the number of protomers that carry the site. Positive values destabilise.

The wild type is prepared once per structure: three Cartesian FastRelax relaxations (`ref2015_cart`, no restraints, five repeats each), of which the lowest-energy one is kept and cached. For each substitution, the residues within 8 Å of the substituted sites form a shell, and a wild-type arm and a mutant arm are each optimised over that shell with `CartesianddG.optimize_structure`.

`--mode single` mutates one protomer only, as a control. With `CB_PER_RESIDUE=1` the script also writes the energy change of each shell residue on each protomer.

## Installation

HC-ddG needs Python 3 and PyRosetta. Apart from PyRosetta it uses only the standard library. It was tested on Ubuntu Linux with Python 3.11 and PyRosetta4.Release.python311.ubuntu.cxx11thread.serialization 2026.21+release.36f0c659eb. Values depend on the PyRosetta build (see Limitations in [details/ADDITIONAL_DETAILS.md](details/ADDITIONAL_DETAILS.md)).

PyRosetta is installed separately under its own licence (see License). One route is the official installer:

    pip install pyrosetta-installer
    python -c "import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()"

Download the repository and run the script from where you like. `homomer_cartesian_ddg.py` is the script. `homomer_cartesian_ddg_annotated.py` is the same code with a short explanation at the top and a one-line comment on each section.

## Usage

Quick check:

    cd example
    bash run_example.sh

This scores three substitutions in the ROP dimer (PDB 1ROP), mutating both protomers, and compares the result with `expected_output.csv`. The first run prepares and caches the relaxed wild type, which can take ten minutes or more. It ends with "Output matches expected_output.csv." on the tested build. A 2026.03 build under Python 3.10 gave 4.18 instead of 3.08 for L41V.

To run your own substitutions:

    python homomer_cartesian_ddg.py structures mutations.csv

`structures` is the folder of structure files and `mutations.csv` is the substitution CSV. The results are written to `hcddg_out.csv`. The input format, output columns, options, advanced settings, differences from the standard protocols and limitations are described in [details/ADDITIONAL_DETAILS.md](details/ADDITIONAL_DETAILS.md).

## Benchmark

Preliminary. The per-protein and substitution-level results are withheld until the accompanying paper is published. Only summary results are given here.

We compared three treatments of a homo-oligomer on the same substitutions: the substitution in one protomer, in every protomer with no symmetry imposed, and in every protomer with symmetry enforced during relaxation. The panel is 20 homo-oligomers with 420 substitutions of measured folding ddG. 1SCE, a four-chain assembly, was excluded from the panel because every method, including ThermoMPNN, failed to rank its substitutions, which points to a problem with its structure or experimental data rather than with any one method. It will be re-examined in the full benchmark. The statistic is the Spearman correlation between predicted and measured ddG within each protein, averaged over proteins, because the practical question is how to rank candidate substitutions inside one protein.

Every program ran on the same structures, including the two trimers and the tetramer. ThermoMPNN was given the whole assembly as structural context, scored the substitution on each protomer separately and the values were averaged over protomers, so it does not model the simultaneous substitution in every chain. Its double-mutant model (ThermoMPNN-D) was not run.

Mean per-protein Spearman correlation with the experimental ddG (preliminary):

| Method | Mean |
|---|---|
| Every protomer (HC-ddG) | 0.531 |
| One protomer | 0.408 |
| Enforced symmetry (19 of the 20 proteins) | 0.187 |
| Faithful reimplementation of Frenz et al. | 0.520 |
| ThermoMPNN | 0.454 |

Preliminary findings:

- Mutating every protomer ranked substitutions better than mutating one.
- Enforcing symmetry cost the most accuracy.
- HC-ddG performed comparably to a faithful reimplementation of the Frenz et al. protocol, which also mutated every chain, and scored above ThermoMPNN.
- Mutating every protomer improved ranking both near and distant from the interface.

The full analysis will be released with the paper.

## License

HC-ddG is released under the Creative Commons Attribution 4.0 International licence (CC BY 4.0; see LICENSE). You may share and adapt it, including commercially, if you give appropriate credit, link to the licence and indicate any changes. To give credit, cite the Zenodo record below.

HC-ddG calls the Rosetta software suite through PyRosetta and contains none of its code. The author claims no rights in Rosetta, PyRosetta or the Cartesian ddG method of Park et al. and Frenz et al., which remain the property of their developers and licensors. Running HC-ddG needs PyRosetta under its own licence, which is free for non-commercial use. Commercial use of PyRosetta needs a licence from the University of Washington. HC-ddG is independent work and is not affiliated with or endorsed by Rosetta Commons or the University of Washington.

`example/structures/1ROP.pdb` holds the coordinates (ATOM records only) of Protein Data Bank entry 1ROP, which the wwPDB releases under CC0 1.0. The structure is from Banner DW, Kokkinidis M, Tsernoglou D. Structure of the ColE1 rop protein at 1.7 Å resolution. J Mol Biol. 1987;196:657-675. https://doi.org/10.1016/0022-2836(87)90039-8

## Citation

Peterson BE, Tsai J, Schacht P. HC-ddG: Homomer Cartesian ddG for homo-oligomer stability (version 1.0.0) [software]. Zenodo. 2026. https://doi.org/10.5281/zenodo.20808718

Please also cite the work HC-ddG builds on:

- Park H, Bradley P, Greisen P Jr, Liu Y, Mulligan VK, Kim DE, Baker D, DiMaio F. Simultaneous optimization of biomolecular energy functions on features from small molecules and macromolecules. J Chem Theory Comput. 2016;12:6201-6212. https://doi.org/10.1021/acs.jctc.6b00819
- Frenz B, Lewis SM, King I, DiMaio F, Park H, Song Y. Prediction of protein mutational free energy: benchmark and sampling improvements increase classification accuracy. Front Bioeng Biotechnol. 2020;8:558247. https://doi.org/10.3389/fbioe.2020.558247
- Kellogg EH, Leaver-Fay A, Baker D. Role of conformational sampling in computing mutation-induced changes in protein structure and stability. Proteins. 2011;79:830-838. https://doi.org/10.1002/prot.22921
- Chaudhury S, Lyskov S, Gray JJ. PyRosetta: a script-based interface for implementing molecular modeling algorithms using Rosetta. Bioinformatics. 2010;26:689-691. https://doi.org/10.1093/bioinformatics/btq007
- Leaver-Fay A, Tyka M, Lewis SM, Lange OF, et al. ROSETTA3: an object-oriented software suite for the simulation and design of macromolecules. Methods Enzymol. 2011;487:545-574. https://doi.org/10.1016/B978-0-12-381270-4.00019-6
- Alford R, Leaver-Fay A, Jeliazkov J, O'Meara M, et al. The Rosetta all-atom energy function for macromolecular modeling and design. J Chem Theory Comput. 2017;13:3031-3048. https://doi.org/10.1021/acs.jctc.7b00125
- Koehler Leman J, Weitzner BD, Lewis SM, Adolf-Bryfogle J, et al. Macromolecular modeling and design in Rosetta: recent methods and frameworks. Nat Methods. 2020;17:665-680. https://doi.org/10.1038/s41592-020-0848-2

Program compared in the benchmark:

- Dieckhaus H, Brocidiacono M, Randolph NZ, Kuhlman B. Transfer learning to leverage larger datasets for improved prediction of protein stability changes. Proc Natl Acad Sci U S A. 2024;121:e2314853121. https://doi.org/10.1073/pnas.2314853121 (ThermoMPNN)
