# OpenMM simulation example


**This branch is an updated version of Tim's original tutorial, with a more recent toolkit release.**
**Only updated-simulate-complex.py has been updated; other Python scripts are from the original repository.**

@Ialibay also contributed to updating this script.

Create an environment:

```
conda create env -f environment.yaml
conda activate example-simulate-complex
```

## updated-simulate-complex.py: Basic Complex Simulation

How to set up a protein-ligand complex simulation.

This is set up with the same defaults as Tim's original script.

```
python updated-simulate-complex.py --help
usage: Simulate a protein-ligand complex, with or without water
       [-h] [--pdb-file PDB_FILE] [--small-mol-file SMALL_MOL_FILE]
       [--output-name OUTPUT_NAME] [--temperature TEMPERATURE] [--solvate]
       [--box-padding BOX_PADDING]
       [--solvent-ionic-strength SOLVENT_IONIC_STRENGTH]
       [--n-equilibration-steps N_EQUILIBRATION_STEPS]
       [--reporting-interval REPORTING_INTERVAL]
       [--n-simulation-steps N_SIMULATION_STEPS]
       [--protein-forcefield PROTEIN_FORCEFIELD]
       [--ligand-forcefield {gaff-2.11,openff-2.0.0,openff-1.3.1,openff-1.3.0}]
       [--solvent-forcefield SOLVENT_FORCEFIELD] [--water-model WATER_MODEL]

optional arguments:
  -h, --help            show this help message and exit
  --pdb-file PDB_FILE   PDB file of protein
  --small-mol-file SMALL_MOL_FILE
                        SDF file of small molecules
  --output-name OUTPUT_NAME
                        Output file name prefix
  --temperature TEMPERATURE
                        Simulation temperature in Kelvin
  --solvate             Whether to solvate the system
  --box-padding BOX_PADDING
                        Padding around protein in nanometers to solvate. This
                        option is ignored if --solvate is not set
  --solvent-ionic-strength SOLVENT_IONIC_STRENGTH
                        Ionic strength of solvent in molar. This option is
                        ignored if --solvate is not set
  --n-equilibration-steps N_EQUILIBRATION_STEPS
                        Number of equilibration steps
  --reporting-interval REPORTING_INTERVAL
                        Reporting interval
  --n-simulation-steps N_SIMULATION_STEPS
                        Number of simulation steps
  --protein-forcefield PROTEIN_FORCEFIELD
                        Protein force field (openmmforcefields)
  --ligand-forcefield {gaff-2.11,openff-2.0.0,openff-1.3.1,openff-1.3.0}
                        Ligand force field
  --solvent-forcefield SOLVENT_FORCEFIELD
                        Solvent force field (openmmforcefields)
  --water-model WATER_MODEL
                        Water model to use for solvation
```

The small molecule file is read using RDKit and then processed to:
* Add hydrogens
* Define the stereochemistry

It can include multiple molecules if the system has cofactors.

The protein is then read in PDB format and added to a Modeller object.
The ligand is then added to the Modeller to generate the complex.
The system is then prepared using the appropriate force fields.

An example is given in `examples/eg5`.

