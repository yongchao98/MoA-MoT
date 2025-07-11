try:
    # This script uses the RDKit library for chemical informatics.
    # If you don't have it installed, you can install it via pip:
    # pip install rdkit-pypi
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
except ImportError:
    print("This script requires the RDKit library.")
    print("Please install it using: pip install rdkit-pypi")
    # In case of an import error, we exit to avoid further errors.
    exit()

# The structure of the final product, Compound A, is (3-hydroxypyridin-2-yl)(phenylamino)acetonitrile.
# We can represent this molecule using a SMILES string.
smiles_A = "Oc1ncccc1C(C#N)Nc1ccccc1"

# Create a molecule object from the SMILES string
mol_A = Chem.MolFromSmiles(smiles_A)

if mol_A is not None:
    print("Analysis of Compound A: (3-hydroxypyridin-2-yl)(phenylamino)acetonitrile")
    
    # Calculate and display the molecular formula and exact molecular weight.
    molecular_formula = CalcMolFormula(mol_A)
    molecular_weight = Descriptors.ExactMolWt(mol_A)
    
    print(f"\nMolecular Formula: {molecular_formula}")
    print(f"Exact Molecular Weight: {molecular_weight:.4f} Da")

    # To fulfill the requirement "output each number in the final equation",
    # we will display the count for each element in the molecular formula.
    print("\nElemental composition of Compound A:")
    
    atom_counts = {}
    for atom in mol_A.GetAtoms():
        symbol = atom.GetSymbol()
        atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
        
    # Print atom counts in the conventional order (C, H, then alphabetically).
    if 'C' in atom_counts:
        print(f"Carbon (C): {atom_counts['C']}")
    if 'H' in atom_counts:
        print(f"Hydrogen (H): {atom_counts['H']}")
    
    # Handle other elements alphabetically.
    other_elements = sorted([symbol for symbol in atom_counts if symbol not in ['C', 'H']])
    for symbol in other_elements:
        element_name = Chem.Atom(symbol).GetPeriodicEntry().name
        print(f"{element_name} ({symbol}): {atom_counts[symbol]}")
        
else:
    print("Error: Could not parse the SMILES string for Compound A.")
