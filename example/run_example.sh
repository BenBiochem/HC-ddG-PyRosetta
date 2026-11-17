#!/bin/bash
# Scores three substitutions in the ROP dimer (PDB 1ROP), mutating both
# protomers, then compares the result with expected_output.csv. Values
# depend on the PyRosetta build; see the README.
cd "$(dirname "$0")"
CB_EXPECT_N=2 CB_WTDIR=wt_cache python ../homomer_cartesian_ddg.py \
	structures mutations.csv --out example_out.csv --workers 3
diff example_out.csv expected_output.csv && echo "Output matches expected_output.csv."
