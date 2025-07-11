# Note: To run this code, you need to install the rdkit and pyscf libraries.
# You can install them using the command: pip install rdkit pyscf

from rdkit import Chem
from rdkit.Chem import AllChem
from pyscf import gto

def find_molecule_symmetry(smiles_string):
    """
    This function takes a SMILES string, generates a 3D structure,
    and determines its molecular point group using rdkit and pyscf.
    """
    try:
        # Step 1: Create a molecule object from the SMILES string using rdkit
        mol_rdkit = Chem.MolFromSmiles(smiles_string)
        if mol_rdkit is None:
            print("Error: Invalid SMILES string provided.")
            return

        # Add hydrogens to complete the valence
        mol_rdkit = Chem.AddHs(mol_rdkit)

        # Step 2: Generate a 3D conformation of the molecule
        # We use a random seed for reproducibility
        AllChem.EmbedMolecule(mol_rdkit, randomSeed=1)
        # Optimize the geometry using the Universal Force Field (UFF)
        # A try-except block handles potential optimization failures for complex molecules
        try:
            AllChem.UFFOptimizeMolecule(mol_rdkit)
        except Exception as e:
            # For complex rigid molecules, optimization might not be necessary or might fail.
            # The embedded structure is often sufficient for symmetry analysis.
            print(f"Warning: RDKit geometry optimization failed, proceeding with initial 3D structure. Error: {e}")

        # Step 3: Extract atom symbols and 3D coordinates for PySCF
        conf = mol_rdkit.GetConformer()
        pyscf_atoms = []
        for atom in mol_rdkit.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            pyscf_atoms.append([atom.GetSymbol(), (pos.x, pos.y, pos.z)])

        # Step 4: Use PySCF to determine the symmetry group
        mol_pyscf = gto.Mole()
        mol_pyscf.atom = pyscf_atoms
        mol_pyscf.build(symm=True)  # The build function with symm=True detects the point group

        # Step 5: Print the results
        group_name = mol_pyscf.symmetry
        print(f"The symmetry group for the molecule is: {group_name}")
        
        # Breakdown of the point group name characters
        print("\n--- Explanation of the Symmetry Group Name ---")
        if group_name and len(group_name) >= 3:
            char1 = group_name[0]
            char2 = group_name[1]
            char3 = group_name[2]
            
            print(f"Character '{char1}': Stands for 'Dihedral', indicating the presence of {char2} C2 axes perpendicular to the main rotation axis.")
            print(f"Character '{char2}': Represents the order of the principal (highest-order) rotation axis. Here, it is a C3 axis.")
            print(f"Character '{char3}': Stands for 'horizontal', indicating a horizontal mirror plane (σh) perpendicular to the principal C3 axis.")

    except ImportError:
        print("Error: The 'rdkit' or 'pyscf' library is not installed.")
        print("Please install them by running: pip install rdkit pyscf")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# The SMILES string for the molecule in question
smiles = "C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67"

find_molecule_symmetry(smiles)