# This script determines the symmetry point group of a molecule.
#
# The original SMILES string provided by the user is syntactically incorrect.
# Based on the visible pattern of three ethynyl (C#C) groups on an aromatic
# system, we assume the intended molecule is 1,3,5-triethynylbenzene,
# which is a common, highly symmetric molecule fitting the description.
#
# This script requires the RDKit and ASE libraries. You can install them using pip:
# pip install rdkit-pypi ase

import sys

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from ase import Atoms
    from ase.symmetry.pointgroup import get_pointgroup
except ImportError:
    print("This script requires RDKit and ASE. Please install them first:")
    print("pip install rdkit-pypi ase")
    sys.exit(1)

def find_molecule_symmetry(smiles_string, molecule_name):
    """
    Calculates and prints the point group of a molecule from its SMILES string.
    """
    print(f"Analyzing molecule: {molecule_name}")
    print(f"SMILES: {smiles_string}\n")

    # Step 1: Create a molecule object from the SMILES string
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        print(f"Error: Could not parse the SMILES string for {molecule_name}.")
        return

    # Step 2: Add hydrogens to the molecule
    mol = Chem.AddHs(mol)

    # Step 3: Generate a 3D conformation for the molecule
    # We use a fixed random seed for reproducibility of the 3D structure.
    status = AllChem.EmbedMolecule(mol, randomSeed=42)
    if status == -1:
        print("Error: Failed to generate 3D coordinates.")
        return

    # Step 4: Optimize the geometry using a force field for a realistic structure
    try:
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        # Fallback to UFF optimizer if MMFF fails
        try:
            AllChem.UFFOptimizeMolecule(mol)
        except Exception as e:
            print(f"Warning: Force field optimization failed: {e}")
            print("Proceeding with unoptimized geometry, which may affect accuracy.")

    # Step 5: Convert the RDKit molecule to an ASE (Atomic Simulation Environment) Atoms object
    try:
        positions = mol.GetConformer().GetPositions()
        atomic_numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
        ase_atoms = Atoms(numbers=atomic_numbers, positions=positions)
    except Exception as e:
        print(f"Error converting molecule to ASE format: {e}")
        return
        
    # Step 6: Use ASE to determine the point group.
    # The tolerance parameter defines the precision for symmetry detection.
    try:
        pg = get_pointgroup(ase_atoms, tolerance=0.1)
        point_group_name = pg.schonflies
        print(f"The determined Schoenflies symbol for the point group is: {point_group_name}")
    except Exception as e:
        print(f"Could not determine point group: {e}")


if __name__ == "__main__":
    # The SMILES string for 1,3,5-triethynylbenzene.
    corrected_smiles = "C#Cc1cc(C#C)cc(C#C)c1"
    
    print("The provided SMILES string is invalid. We will analyze a likely alternative.\n")
    find_molecule_symmetry(corrected_smiles, "1,3,5-triethynylbenzene")