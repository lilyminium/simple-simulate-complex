#!/usr/bin/env bash

source ~/.bashrc
conda activate example-simulate-complex

python ../../updated-simulate-complex.py                    \
    --pdb-file                  eg5_protein.pdb             \
    --small-mol-file            eg5_small_mols.sdf          \
    --output-name               output                      \
    --n-equilibration-steps     5                           \
    --n-simulation-steps        10                          \
    --solvate                                               \
    --solvent-ionic-strength    0.15                        \
    --box-padding               1.0                         \
    --protein-forcefield        "amber/ff14SB.xml"          \
    --ligand-forcefield         "openff-2.0.0"              