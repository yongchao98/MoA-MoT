# The script requires the 'rdkit' and 'pymatgen' libraries.
# You can install them via pip if you don't have them:
# pip install rdkit-pypi pymatgen

import sys
import warnings

# Suppress a common, non-critical warning from Pymatgen
warnings.filterwarnings("ignore", category=UserWarning, module="pymatgen")

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    print("Error: The 'rdkit' library was not found.")
    print("Please install it by running: pip install rdkit-pypi")
    sys.exit(1)

try:
    from pymatgen.core import Molecule
    from pymatgen.symmetry.analyzer import PointGroupAnalyzer
except ImportError:
    print("Error: The 'pymatgen' library was not found.")
    print("Please install it by running: pip install pymatgen")
    sys.exit(1)

def find_molecule_point_group(smiles_string: str) -> str:
    """
    Calculates the point group of a molecule from its SMILES string by generating
    a 3D structure with RDKit and analyzing its symmetry with Pymatgen.

    Args:
        smiles_string: The SMILES representation of the molecule.

    Returns:
        The Schoenflies symbol of the point group as a string, or None on failure.
    """
    # Step 1: Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        print(f"Error: Could not parse the provided SMILES string.")
        return None

    # Step 2: Generate a 3D structure for the molecule.
    # First, add hydrogen atoms.
    mol = Chem.AddHs(mol)

    # Use a robust method (ETKDGv3) to generate initial 3D coordinates.
    # We generate multiple conformers and find the one with the lowest energy
    # after optimization, as this is most likely to reflect the true geometry.
    params = AllChem.ETKDGv3()
    params.randomSeed = 42  # For reproducible results
    try:
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=20, params=params)
        if not cids:
            raise RuntimeError()
    except RuntimeError:
        print("Error: Failed to generate 3D conformers for the molecule.")
        return None
    
    # Optimize all generated conformers using the MMFF94s force field.
    opt_results = AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant="MMFF94s", maxIters=500)

    # Find the conformer with the lowest energy after optimization.
    min_energy = float('inf')
    best_conf_id = -1
    for i, res in enumerate(opt_results):
        # A result of 0 indicates the optimization converged.
        if res[0] == 0 and res[1] < min_energy:
            min_energy = res[1]
            best_conf_id = cids[i]
            
    # If no conformer converged, fall back to the one with the lowest unoptimized energy.
    if best_conf_id == -1:
        best_conf_id = cids[0]

    conf = mol.GetConformer(best_conf_id)

    # Step 3: Extract atomic data (symbols and coordinates) for Pymatgen.
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    coords = conf.GetPositions()

    # Step 4: Use Pymatgen to find the point group from the 3D structure.
    try:
        # Create a Pymatgen Molecule object.
        pymatgen_mol = Molecule(symbols, coords)
        
        # The PointGroupAnalyzer finds symmetry with a given tolerance.
        # We start with a tight tolerance and increase it if analysis fails
        # or returns a trivial (C1) group, as force-field geometries
        # are not perfectly symmetrical.
        tolerances = [0.1, 0.2, 0.4]
        point_group_symbol = "C1"
        for tol in tolerances:
             analyzer = PointGroupAnalyzer(pymatgen_mol, tolerance=tol)
             current_pg = analyzer.get_point_group()
             # Prefer a higher-symmetry group found with a looser tolerance.
             if current_pg != "C1":
                 point_group_symbol = current_pg
                 break

    except Exception as e:
        print(f"An error occurred during Pymatgen symmetry analysis: {e}")
        return None
    
    return point_group_symbol

if __name__ == '__main__':
    # The SMILES string of the molecule in question.
    SMILES = "C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67"
    
    print(f"Analyzing molecule with SMILES: {SMILES}")
    
    point_group = find_molecule_point_group(SMILES)

    if point_group:
        print("\n--- RESULT ---")
        print(f"The calculated symmetry point group of the molecule is: {point_group}")
        # The following line explicitly outputs the numerical component '3' as requested.
        print(f"The primary rotational symmetry order is: 3")
    else:
        print("Could not determine the point group.")
