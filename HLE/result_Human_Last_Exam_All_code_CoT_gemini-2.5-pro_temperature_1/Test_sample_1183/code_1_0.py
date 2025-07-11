import sys
import warnings

# Suppress Pymatgen and RDKit warnings for cleaner output
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning)

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("Error: RDKit is not installed. Please install it using 'pip install rdkit-pypi'", file=sys.stderr)
    sys.exit(1)

try:
    import pymatgen.core as mg
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer
except ImportError:
    print("Error: Pymatgen is not installed. Please install it using 'pip install pymatgen'", file=sys.stderr)
    sys.exit(1)

def find_molecule_point_group(smiles_string: str) -> str:
    """
    Calculates the point group of a molecule from its SMILES string by generating
    a 3D structure with RDKit and analyzing its symmetry with Pymatgen.

    Args:
        smiles_string: The SMILES representation of the molecule.

    Returns:
        The Schoenflies symbol of the point group (e.g., 'C2v', 'D3h').
        Returns an error message if the process fails.
    """
    # 1. Create an RDKit molecule object from the SMILES string.
    mol = Chem.MolFromSmiles(smiles_string)
    if not mol:
        return "Failed: Invalid SMILES string provided."

    # 2. Add hydrogen atoms to the molecule.
    mol = Chem.AddHs(mol)

    # 3. Generate multiple 3D conformers to find a low-energy structure.
    # This increases the chances of finding the most symmetric conformation.
    num_conformers = 50
    p = AllChem.ETKDGv3()
    p.randomSeed = 42 # for reproducibility
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, params=p)
    
    if not cids:
        # Fallback to a single embedding if multiple conformer generation fails
        cid = AllChem.EmbedMolecule(mol, randomSeed=42)
        if cid == -1:
            return "Failed: Could not generate 3D coordinates."

    # 4. Optimize the geometries of all conformers using the MMFF94s force field.
    try:
        opt_results = AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant='MMFF94s', maxIters=2000)
        if not opt_results:
             return "Failed: Geometry optimization failed."
        # Find the index of the conformer with the lowest energy.
        min_energy_idx = min(opt_results, key=lambda x: x[1])[0]
    except Exception:
         # Some complex structures can cause the optimizer to fail.
         min_energy_idx = 0


    # 5. Extract atomic symbols and coordinates from the lowest-energy conformer.
    conformer = mol.GetConformer(min_energy_idx)
    positions = conformer.GetPositions()
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]

    # 6. Create a Pymatgen Molecule object from the symbols and coordinates.
    pmg_molecule = mg.Molecule(symbols, positions)

    # 7. Use Pymatgen's PointGroupAnalyzer to determine the symmetry.
    # A tolerance value is needed to account for small deviations from perfect symmetry
    # in the force field-optimized structure.
    try:
        analyzer = PointGroupAnalyzer(pmg_molecule, tolerance=0.3)
        point_group_symbol = analyzer.get_point_group()
        return point_group_symbol
    except Exception as e:
        return f"Failed: Pymatgen symmetry analysis failed. Error: {e}"

if __name__ == '__main__':
    # The SMILES string provided in the task.
    smiles_input = "C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67"

    # Calculate and print the point group.
    point_group = find_molecule_point_group(smiles_input)
    
    print(point_group)
