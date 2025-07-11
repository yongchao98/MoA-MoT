# To run this code, you may need to install the rdkit library:
# pip install rdkit-pypi

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
except ImportError:
    print("RDKit library not found.")
    print("Please install it using: pip install rdkit-pypi")
    exit()

# Step 1: Identify the product of the reaction.
# The Favorskii rearrangement of 2-chlorocyclohexanone results in a ring contraction.
# The product is cyclopentanecarboxylic acid.
product_name = "cyclopentanecarboxylic acid"
smiles_string = "O=C(O)C1CCCC1"
mol = Chem.MolFromSmiles(smiles_string)

if mol is None:
    print(f"Error: Could not create molecule from SMILES string: {smiles_string}")
else:
    # Step 2: Calculate basic molecular properties that contribute to complexity.
    # The BertzCT index is a sophisticated value derived from the graph theory of the molecule.
    # We list some of the underlying properties here for context.
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    num_atoms = mol.GetNumAtoms()
    num_bonds = mol.GetNumBonds()
    num_rings = Descriptors.RingCount(mol)

    # Step 3: Calculate the Böttcher (Bertz) Molecular Complexity.
    # This is a single, comprehensive index. It's not a simple sum of the parts above.
    complexity_index = Descriptors.BertzCT(mol)

    # Step 4: Display the results.
    print(f"The product of the Favorskii rearrangement of 2-chlorocyclohexanone is {product_name}.")
    print("\nMolecular properties contributing to complexity:")
    print(f"Number of non-hydrogen atoms = {num_heavy_atoms}")
    print(f"Total number of atoms = {num_atoms}")
    print(f"Total number of bonds = {num_bonds}")
    print(f"Number of rings = {num_rings}")
    
    print("\nThe final Böttcher Molecular Complexity is calculated using the BertzCT index.")
    print(f"Final Value = {complexity_index:.4f}")
