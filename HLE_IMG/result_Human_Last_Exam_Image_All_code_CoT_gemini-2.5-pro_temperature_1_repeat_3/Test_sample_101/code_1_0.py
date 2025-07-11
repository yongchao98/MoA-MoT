# First, ensure you have RDKit installed: pip install rdkit-pypi
try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    from collections import Counter
except ImportError:
    print("RDKit library is required for this script.")
    print("Please install it using: pip install rdkit-pypi")
    exit()

# --- Chemical Analysis ---
# Step 1: 3-hydroxy-pyridine-2-carbaldehyde + Aniline -> Imine intermediate + H2O
# Step 2: Imine intermediate + NaCN -> Compound A (alpha-aminonitrile)
# This is a Strecker-type reaction where a cyanide ion adds across the C=N double bond.

# --- Compound A Details ---
# The final product is 2-((3-hydroxypyridin-2-yl)(phenylamino))acetonitrile.
# We can represent it using a SMILES string.
smiles_compound_A = "N#CC(Nc1ccccc1)c1nccc(O)c1"
compound_A = Chem.MolFromSmiles(smiles_compound_A)

print("--- Analysis of the Final Product (Compound A) ---")

if compound_A:
    # RDKit automatically handles implicit hydrogens, but we make them explicit for counting.
    compound_A_with_H = Chem.AddHs(compound_A)

    # Calculate the molecular formula
    molecular_formula = rdMolDescriptors.CalcMolFormula(compound_A_with_H)

    # Count each type of atom to get the 'numbers in the final equation'
    atom_symbols = [atom.GetSymbol() for atom in compound_A_with_H.GetAtoms()]
    atom_counts = Counter(atom_symbols)

    print(f"Product Name: 2-((3-hydroxypyridin-2-yl)(phenylamino))acetonitrile")
    print(f"SMILES String: {smiles_compound_A}")
    print(f"\nMolecular Formula: {molecular_formula}")
    print("\nAtom Counts (Numbers in the molecular formula):")
    print(f"Number of Carbon (C) atoms: {atom_counts['C']}")
    print(f"Number of Hydrogen (H) atoms: {atom_counts['H']}")
    print(f"Number of Nitrogen (N) atoms: {atom_counts['N']}")
    print(f"Number of Oxygen (O) atoms: {atom_counts['O']}")
else:
    print("Error: Could not process the chemical structure from the SMILES string.")
