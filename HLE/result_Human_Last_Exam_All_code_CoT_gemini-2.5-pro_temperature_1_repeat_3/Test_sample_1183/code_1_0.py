import sys
from rdkit import Chem
from rdkit.Chem import AllChem

# Suppress PySCF warning messages for a cleaner output
class SuppressStream:
    def __enter__(self):
        self._original_stdout = sys.stdout
        self._original_stderr = sys.stderr
        sys.stdout = self
        sys.stderr = self

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self._original_stdout
        sys.stderr = self._original_stderr

    def write(self, text):
        pass

    def flush(self):
        pass

try:
    with SuppressStream():
        from pyscf import gto
        from pyscf.symm import addons
except ImportError:
    print("Error: PySCF or RDKit is not installed.")
    print("Please install them by running: pip install pyscf rdkit-pypi")
    sys.exit(1)


def find_point_group(smiles_string):
    """
    Calculates the point group of a molecule from its SMILES string.

    Args:
        smiles_string (str): The SMILES representation of the molecule.

    Returns:
        str: The Schoenflies symbol of the point group, or an error message.
    """
    try:
        # 1. Create RDKit molecule object from SMILES and add hydrogens
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            return "Error: Invalid SMILES string."
        mol = Chem.AddHs(mol)

        # 2. Generate 3D coordinates and optimize the geometry
        # Use a fixed random seed for reproducibility
        AllChem.EmbedMolecule(mol, randomSeed=42)
        # Use MMFF94 force field for optimization
        if AllChem.MMFFOptimizeMolecule(mol) == -1:
             # If MMFF fails, try UFF
             AllChem.UFFOptimizeMolecule(mol)


        # 3. Prepare atom data for PySCF
        conformer = mol.GetConformer()
        pyscf_atoms = []
        for atom in mol.GetAtoms():
            pos = conformer.GetAtomPosition(atom.GetIdx())
            pyscf_atoms.append([atom.GetSymbol(), (pos.x, pos.y, pos.z)])

        # 4. Create PySCF molecule object and detect symmetry
        # A tolerance of 0.1 Angstrom is a reasonable default
        # for geometries optimized with force fields.
        with SuppressStream():
             mol_pyscf = gto.Mole()
             mol_pyscf.atom = pyscf_atoms
             mol_pyscf.build(dump_input=False, verbose=0)
             point_group = addons.detect_symm(mol_pyscf.atom, tol=0.1)
        
        return point_group

    except Exception as e:
        return f"An error occurred: {e}"

# The SMILES string provided by the user
smiles = 'C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67'

# Get and print the point group
symmetry_group = find_point_group(smiles)
print(symmetry_group)