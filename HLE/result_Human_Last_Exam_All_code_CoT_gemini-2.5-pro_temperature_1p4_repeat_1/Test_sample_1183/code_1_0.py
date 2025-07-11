# In order to run this code, you need to install the rdkit and ase libraries.
# You can install them using pip:
# pip install rdkit ase

import ase
from ase.symmetry import get_point_group
from rdkit import Chem
from rdkit.Chem import AllChem

def find_molecule_symmetry(smiles_string):
    """
    This function takes a SMILES string, generates a 3D structure,
    and determines its molecular point group using rdkit and ase.
    """
    try:
        # Step 1 & 2: Parse SMILES and generate a 3D structure with rdkit
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            print(f"Error: Could not parse SMILES string: {smiles_string}")
            return

        # Add hydrogens
        mol = Chem.AddHs(mol)

        # Generate 3D coordinates using the ETKDG method
        # A random seed is used for reproducibility
        AllChem.EmbedMolecule(mol, AllChem.ETKDG(), randomSeed=1)

        # Optimize the geometry using the Universal Force Field (UFF)
        AllChem.UFFOptimizeMolecule(mol)

        # Step 3: Convert the rdkit molecule to an ase Atoms object
        conformer = mol.GetConformer()
        positions = conformer.GetPositions()
        atomic_numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
        
        atoms = ase.Atoms(numbers=atomic_numbers, positions=positions)
        atoms.center() # Center the molecule at the origin

        # Step 4: Use ase to find the point group
        pg = get_point_group(atoms, tolerance=0.1)
        
        # Step 5: Print the result
        schoenflies_symbol = pg.sch_symbol
        print(f"The symmetry group of the molecule is: {schoenflies_symbol}")

    except ImportError:
        print("Error: RDKit or ASE library not found.")
        print("Please install them using: pip install rdkit ase")
    except Exception as e:
        print(f"An error occurred: {e}")

# The SMILES string for the molecule in question
smiles = "C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67"

# Run the analysis
find_molecule_symmetry(smiles)
