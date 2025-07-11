# The following script requires the RDKit library.
# You can install it using pip:
# pip install rdkit-pypi

from rdkit import Chem
from collections import Counter

def get_atom_counts(mol_smiles):
    """Creates an RDKit molecule, adds hydrogens, and returns a dictionary of atom counts."""
    mol = Chem.MolFromSmiles(mol_smiles)
    if mol:
        mol = Chem.AddHs(mol)  # Add explicit hydrogens to the graph to count them
        return Counter(atom.GetSymbol() for atom in mol.GetAtoms())
    return Counter()

def format_detailed_formula(counts):
    """Formats the atom count dictionary into a detailed 'Symbol(Count)' string."""
    # Define standard atom order for chemical formulas (C, H, then alphabetical)
    order = ['C', 'H'] + sorted([symbol for symbol in counts if symbol not in ['C', 'H']])
    parts = []
    for symbol in order:
        if symbol in counts:
            count = counts[symbol]
            parts.append(f"{symbol}({count})")
    return " ".join(parts)

# 1. Define the reactants and the known product via SMILES strings.
butadiene_smiles = 'C=CC=C'
dienophile_smiles = 'F/C(F)=C(Cl)Cl'
# The product is 4,4-dichloro-3,3-difluorocyclohex-1-ene. Its SMILES is:
product_smiles = 'C1=CC(F)(F)C(Cl)(Cl)CC1'

# 2. Get the atom counts for each molecule.
butadiene_counts = get_atom_counts(butadiene_smiles)
dienophile_counts = get_atom_counts(dienophile_smiles)
product_counts = get_atom_counts(product_smiles)

# 3. Format the detailed formulas.
butadiene_formula_str = format_detailed_formula(butadiene_counts)
dienophile_formula_str = format_detailed_formula(dienophile_counts)
product_formula_str = format_detailed_formula(product_counts)

# The reaction coefficients are all 1.
coeff_reactant1 = 1
coeff_reactant2 = 1
coeff_product = 1

# 4. Print the final result.
print("The product of the reaction between butadiene and 1,1-dichloro-2,2-difluoroethene is:")
print("4,4-dichloro-3,3-difluorocyclohex-1-ene")
print("\nThe balanced chemical equation with a breakdown of each atom count is:")
print(
    f"({coeff_reactant1}) {butadiene_formula_str}   +   ({coeff_reactant2}) {dienophile_formula_str}   ->   ({coeff_product}) {product_formula_str}"
)
