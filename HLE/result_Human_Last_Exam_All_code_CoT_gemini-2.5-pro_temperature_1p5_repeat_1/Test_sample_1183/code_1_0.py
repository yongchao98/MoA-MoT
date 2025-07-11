import sys
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from ase import Atoms
    from ase.symmetry.analyzer import PointGroupAnalyzer
except ImportError:
    print("This script requires the RDKit and ASE libraries.")
    print("Please install them using pip: pip install rdkit ase")
    sys.exit(1)

def find_molecule_point_group(smiles_string: str):
    """
    Calculates the point group of a molecule from its SMILES string.

    Args:
        smiles_string: The SMILES representation of the molecule.
    """
    print(f"Analyzing molecule with SMILES string: {smiles_string}")

    # Step 1: Create a molecule object from SMILES using RDKit
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        print("Error: Invalid SMILES string provided.")
        return

    # Add hydrogens, as they are part of the molecule's structure and symmetry
    mol = Chem.AddHs(mol)

    # Step 2: Generate a 3D conformation and optimize it
    # Use a fixed random seed for reproducibility
    if AllChem.EmbedMolecule(mol, randomSeed=42) == -1:
        print("Error: Could not generate a 3D conformation.")
        return
    
    try:
        # MMFF94 is a good general-purpose force field for organic molecules
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception as e:
        # The optimization might fail for complex or strained molecules
        print(f"Warning: Geometry optimization failed: {e}. Using unoptimized geometry.")

    # Step 3: Convert the RDKit molecule to an ASE Atoms object
    try:
        # Extract atomic numbers and coordinates from the RDKit conformer
        atomic_numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
        conformer = mol.GetConformer(0)
        positions = conformer.GetPositions()

        # Create the ASE Atoms object
        atoms_obj = Atoms(numbers=atomic_numbers, positions=positions)
    except Exception as e:
        print(f"Error during conversion from RDKit to ASE object: {e}")
        return

    # Step 4: Use ASE to determine the point group
    try:
        # A tolerance is needed to account for small numerical inaccuracies
        # in the generated 3D coordinates. 0.1 Angstrom is a reasonable value.
        tolerance = 0.1
        analyzer = PointGroupAnalyzer(atoms_obj, tolerance=tolerance)
        point_group = analyzer.get_pointgroup()
        
        # Step 5: Print the result
        print("\n--- Result ---")
        print(f"The calculated point group of the molecule is: {point_group}")

    except Exception as e:
        print(f"Error during symmetry analysis: {e}")


if __name__ == "__main__":
    # The SMILES string of the molecule in question
    smiles = "C#Cc1cc2ccc3c(C#C)cc4ccc5c(C#C)cc6ccc1c7c2c3c4c5c67"
    find_molecule_point_group(smiles)