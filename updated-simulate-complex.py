#!/usr/bin/env python3

import argparse
import logging
import sys
import time
import typing

import openmm
import openmm.app
from openff.toolkit.topology import Molecule
from openff.units import unit
from openff.units.openmm import to_openmm

if typing.TYPE_CHECKING:
    import openmm

logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(
    "Simulate a protein-ligand complex, with or without water"
)
parser.add_argument(
    "--pdb-file", type=str, default="protein.pdb", help="PDB file of protein"
)
parser.add_argument(
    "--small-mol-file", type=str, default="ligand.sdf", help="SDF file of small molecules"
)
parser.add_argument(
    "--output-name", type=str, default="output", help="Output file name prefix"
)
parser.add_argument(
    "--temperature", type=int, default=300, help="Simulation temperature in Kelvin"
)
parser.add_argument("--solvate", action="store_true", help="Whether to solvate the system")
parser.add_argument(
    "--box-padding",
    type=float,
    default=1.0,
    help=(
        "Padding around protein in nanometers to solvate. "
        "This option is ignored if --solvate is not set"
    ),
)
parser.add_argument(
    "--solvent-ionic-strength",
    type=float,
    default=0.15,
    help=(
        "Ionic strength of solvent in molar. "
        "This option is ignored if --solvate is not set"
    ),
)
parser.add_argument(
    "--n-equilibration-steps",
    type=int,
    default=20,
    help="Number of equilibration steps",
)
parser.add_argument(
    "--reporting-interval", type=int, default=1000, help="Reporting interval"
)
parser.add_argument(
    "--n-simulation-steps", type=int, default=500, help="Number of simulation steps"
)
parser.add_argument(
    "--protein-forcefield",
    type=str,
    default="amber/ff14SB.xml",
    help="Protein force field (openmmforcefields)",
)
parser.add_argument(
    "--ligand-forcefield",
    type=str,
    choices=["gaff-2.11", "openff-2.0.0", "openff-1.3.1", "openff-1.3.0"],
    default="gaff-2.11",
    help="Ligand force field",
)
parser.add_argument(
    "--solvent-forcefield",
    type=str,
    default="amber/tip3p_standard.xml",
    help="Solvent force field (openmmforcefields)",
)
parser.add_argument(
    "--water-model",
    type=str,
    default="tip3p",
    help="Water model to use for solvation",
)

def prepare_ligand(ligand_file: str) -> Molecule:
    """
    Read the molfile into RDKit, add Hs, and assign stereochemistry.
    Returns an OpenFF Molecule object

    Parameters
    ----------
    ligand_file : str
        Path to ligand SDF file

    Returns
    -------
    ligand_mols : list[openff.toolkit.topology.Molecule]
    """
    from rdkit import Chem
    supp = Chem.SDMolSupplier(ligand_file, removeHs=False)
    ligand_mols = []
    for lig in supp:
        offmol = Molecule.from_rdkit(lig, allow_undefined_stereo=True)
        ligand_mols.append(offmol)
        logger.info(f"Prepared ligand with {offmol.n_atoms} atoms")
    return ligand_mols


def create_openmm_system_with_openmm_forcefields(
    pdb_file: str,
    ligand_molecules: list[Molecule],
    solvate: bool = True,
    box_padding: float = 1.0,  # nanometers
    solvent_ionic_strength: float = 0.15,  # molar
    output_name: str = "output",
    protein_forcefield: str = "amber/ff14SB.xml",
    ligand_forcefield: str = "gaff-2.11",
    solvent_forcefield: str = "amber/tip3p_standard.xml",
    water_model: str = "tip3p",
) -> "openmm.System":
    """
    Create OpenMM system using only OpenMM force fields

    Parameters
    ----------
    pdb_file : str
        Path to protein PDB file to use. This may or may
        not be solvated.
    ligand_molecules : list(openff.toolkit.topology.Molecule)
        List of ligand OpenFF Molecules
    solvate : bool, optional
        Whether or not to solvate the complex
    box_padding : float, optional
        Padding around protein in nanometers to solvate.
        This option is ignored if ``solvate`` is False.
    solvent_ionic_strength : float, optional
        Ionic strength of solvent in molar.
        This option is ignored if ``solvate`` is False.
    output_name : str, optional
        Prefix for output files, by default "output"
    protein_forcefield : str, optional
        Protein force field (openmmforcefields), by default "amber/ff14SB.xml"
    ligand_forcefield : str, optional
        Ligand force field, by default "gaff-2.11"
    solvent_forcefield : str, optional
        Solvent force field (openmmforcefields), by default "amber/tip3p_standard.xml"
    water_model : str, optional
        Water model, by default "tip3p" -- take care this matches the force field.

    Returns
    -------
    openmm_system : openmm.System
        OpenMM system of complex, solvated or not
    modeller : openmm.app.Modeller
    """
    from openmmforcefields.generators import SystemGenerator

    forcefield_kwargs = dict(
        constraints=openmm.app.HBonds,
        rigidWater=True,
        removeCMMotion=False,
        hydrogenMass=3 * openmm.unit.amu,
    )

    generator = SystemGenerator(
        forcefields=[protein_forcefield, solvent_forcefield],
        small_molecule_forcefield=ligand_forcefield,
        forcefield_kwargs=forcefield_kwargs,
        molecules=ligand_molecules
    )

    protein_pdb = openmm.app.PDBFile(pdb_file)
    logger.info(f"Found protein {len(protein_pdb.positions)} atoms")

    modeller = openmm.app.Modeller(protein_pdb.topology, protein_pdb.positions)
    for mol in ligand_molecules:
        modeller.add(
            mol.to_topology().to_openmm(),
            to_openmm(mol.conformers[0]),
        )

    output_file = f"{output_name}_complex.pdb"
    with open(output_file, "w") as f:
        openmm.app.PDBFile.writeFile(
            modeller.topology,
            modeller.positions,
            f,
        )
    logger.info(f"Created complex with {len(modeller.positions)} atoms")
    logger.info(f"Saved complex to {output_file}")

    if solvate:
        modeller.addSolvent(
            generator.forcefield,
            model=water_model,
            ionicStrength=solvent_ionic_strength * openmm.unit.molar,
            padding=box_padding * openmm.unit.nanometer,
        )
        logger.info(f"Solvated system for total {len(modeller.positions)} atoms")

        output_file = f"{output_name}_complex_solvated.pdb"
        with open(output_file, "w") as f:
            openmm.app.PDBFile.writeFile(
                modeller.topology,
                modeller.positions,
                f,
            )
        logger.info(f"Saved solvated complex to {output_file}")

    system = generator.create_system(
        modeller.topology,
        molecules=ligand_molecules,
    )
    logger.info(f"OpenMM system has {len(modeller.positions)} atoms")
    return system, modeller


def set_precision() -> "openmm.Platform":
    """
    Set simulation precision and get quickest platform

    Returns
    -------
    platform : openmm.Platform
        The platform to use for simulation
    """
    current_speed = 0
    for i in range(openmm.Platform.getNumPlatforms()):
        p = openmm.Platform.getPlatform(i)
        if p.getSpeed() > current_speed:
            current_speed = p.getSpeed()
            platform = p

    if platform.getName() == "CUDA" or platform.getName() == "OpenCL":
        platform.setPropertyDefaultValue("Precision", "mixed")
        logger.info(f"Set precision for platform {platform.getName()} to mixed")
    return platform


def simulate(
    modeller: "openmm.app.Modeller",
    openmm_system: "openmm.System",
    platform: "openmm.Platform",
    temperature: float = 300,  # kelvin
    output_name: str = "output",
    n_equilibration_steps: int = 200,
    reporting_interval: int = 1000,
    n_simulation_steps: int = 5000,
):
    """
    Setup and run simulation with OpenMM


    Parameters
    ----------
    openff_topology : openff.toolkit.topology.Topology
        Topology of complex, solvated or not
    openmm_system : openmm.System
        OpenMM system of complex, solvated or not
    platform : openmm.Platform
        OpenMM platform to use
    temperature : float, optional
        Simulation temperature in Kelvin, by default 300
    output_name : str, optional
        Prefix for output files, by default "output"
    n_equilibration_steps : int, optional
        Number of equilibration steps, by default 200
    reporting_interval : int, optional
        Reporting interval, by default 1000
    n_simulation_steps : int, optional
        Number of simulation steps, by default 5000

    """
    integrator = openmm.LangevinIntegrator(
        temperature * openmm.unit.kelvin,
        1.0 / openmm.unit.picoseconds,
        0.002 * openmm.unit.picoseconds,
    )
    logger.info(f"Uses periodic box: {openmm_system.usesPeriodicBoundaryConditions()}")
    logger.info(f"Default periodic box: {openmm_system.getDefaultPeriodicBoxVectors()}")

    simulation = openmm.app.Simulation(
        modeller.topology,
        openmm_system,
        integrator,
        platform=platform,
    )
    simulation.context.setPositions(modeller.positions)
    logger.info(
        f"Starting energy: {simulation.context.getState(getEnergy=True).getPotentialEnergy()}"
    )

    logger.info("Minimizing")
    simulation.minimizeEnergy()

    min_pdb_file = output_name + "_minimized.pdb"
    positions = simulation.context.getState(
        getPositions=True, enforcePeriodicBox=False
    ).getPositions()
    with open(min_pdb_file, "w") as f:
        openmm.app.PDBFile.writeFile(
            simulation.topology, positions, file=f, keepIds=True
        )
    logger.info(f"Saved minimized structure to {min_pdb_file}")

    simulation.context.setVelocitiesToTemperature(temperature * openmm.unit.kelvin)
    logger.info("Equilibrating")
    simulation.step(n_equilibration_steps)

    output_traj_dcd = output_name + "_traj.dcd"
    # output_traj_pdb = output_name + "_traj.pdb"
    simulation.reporters.append(
        openmm.app.DCDReporter(
            output_traj_dcd, reporting_interval, enforcePeriodicBox=False
        )
    )
    simulation.reporters.append(
        openmm.app.StateDataReporter(
            sys.stdout,
            reporting_interval * 5,
            step=True,
            potentialEnergy=True,
            temperature=True,
        )
    )
    logger.info(f"Starting simulation with {n_simulation_steps} steps")
    t0 = time.time()
    simulation.step(n_simulation_steps)
    t1 = time.time()
    logger.info(f"Completed simulation in {t1 - t0} seconds")


def simulate_complex(
    pdb_file: str = "protein.pdb",
    ligand_file: str = "ligand1.mol",
    output_name: str = "output",
    temperature: int = 300,
    solvate: bool = True,
    box_padding: float = 1.0,  # nanometers
    solvent_ionic_strength: float = 0.15,  # molar
    n_equilibration_steps: int = 20,
    reporting_interval: int = 1000,
    n_simulation_steps: int = 500,
    protein_forcefield: str = "amber/ff14SB.xml",
    ligand_forcefield: str = "gaff-2.11",
    solvent_forcefield: str = "amber/tip3p_standard.xml",
    water_model: str = "tip3p",
):
    """
    Simulate a protein-ligand complex, with and without solvation

    Parameters
    ----------
    pdb_file : str, optional
        Path to protein PDB file, by default "protein.pdb"
    ligand_file : str, optional
        Path to ligand MOL file, by default "ligand1.mol"
    output_name : str, optional
        Prefix for output files, by default "output"
    temperature : int, optional
        Simulation temperature in Kelvin, by default 300
    solvate : bool, optional
        Whether or not to solvate the protein
    box_padding : float, optional
        Padding around protein in nanometers to solvate.
        This option is ignored if ``solvate`` is False.
    solvent_ionic_strength : float, optional
        Ionic strength of solvent in molar.
        This option is ignored if ``solvate`` is False.
    n_equilibration_steps : int, optional
        Number of equilibration steps, by default 20
    reporting_interval : int, optional
        Reporting interval, by default 1000
    n_simulation_steps : int, optional
        Number of simulation steps, by default 500
    protein_forcefield : str, optional
        Protein force field (openmmforcefields), by default "amber/ff14SB.xml"
    ligand_forcefield : str, optional
        Ligand force field, by default "gaff-2.11"
    solvent_forcefield : str, optional
        Solvent force field (openmmforcefields), by default "amber/tip3p_standard.xml"
    water_model : str, optional
        Water model to use for solvation, by default "tip3p".
        Take care that this matches the force field.
    """
    ligands = prepare_ligand(ligand_file)
    system, modeller = create_openmm_system_with_openmm_forcefields(
        pdb_file,
        ligands,
        solvate,
        box_padding=box_padding,
        solvent_ionic_strength=solvent_ionic_strength,
        protein_forcefield=protein_forcefield,
        ligand_forcefield=ligand_forcefield,
        solvent_forcefield=solvent_forcefield,
        water_model=water_model,
    )

    simulate(
        modeller,
        system,
        platform=set_precision(),
        temperature=temperature,
        output_name=output_name,
        n_equilibration_steps=n_equilibration_steps,
        reporting_interval=reporting_interval,
        n_simulation_steps=n_simulation_steps,
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    args = parser.parse_args()
    simulate_complex(
        pdb_file=args.pdb_file,
        ligand_file=args.small_mol_file,
        output_name=args.output_name,
        temperature=args.temperature,
        solvate=args.solvate,
        box_padding=args.box_padding,
        solvent_ionic_strength=args.solvent_ionic_strength,
        n_equilibration_steps=args.n_equilibration_steps,
        reporting_interval=args.reporting_interval,
        n_simulation_steps=args.n_simulation_steps,
        protein_forcefield=args.protein_forcefield,
        ligand_forcefield=args.ligand_forcefield,
        solvent_forcefield=args.solvent_forcefield,
        water_model=args.water_model,
    )
