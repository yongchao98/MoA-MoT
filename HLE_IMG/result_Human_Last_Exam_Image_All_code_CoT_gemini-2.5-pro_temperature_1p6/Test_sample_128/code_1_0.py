import sys

# This script identifies a chemical compound formed in a reaction
# and calculates its key properties using the RDKit library.

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    from rdkit.Chem.rdchem import PeriodicTable
except ImportError:
    print("This script requires the RDKit library to run.")
    print("Please install it, for example using: pip install rdkit-pypi")
    sys.exit(1)

# The reaction shown is a three-component condensation between 2-aminopyridine,
# o-phthalaldehyde, and a cyanide source (TMSCN/KF).
# This reaction is a variant of the Strecker synthesis that results in a
# cyclized product. The product, Compound A, is 2-(pyridin-2-yl)isoindoline-1-carbonitrile.

# We define Compound A using its SMILES (Simplified Molecular-Input Line-Entry System) string.
product_smiles = "N#CC1N(c2ncccc2)Cc2ccccc12"

# Create an RDKit molecule object from the SMILES string
mol = Chem.MolFromSmiles(product_smiles)

if mol is None:
    print(f"Error: Could not process the molecular structure from the SMILES string '{product_smiles}'")
    sys.exit(1)
    
print("--- Analysis of the Final Product (Compound A) ---")
print("Compound A is identified as: 2-(pyridin-2-yl)isoindoline-1-carbonitrile\n")
print(f"SMILES String: {product_smiles}")

# Calculate and display the molecular formula
formula = rdMolDescriptors.CalcMolFormula(mol)
print(f"Molecular Formula: {formula}\n")

print("--- Detailed Composition and Mass Calculation ---")
# Get the periodic table object to access element properties
pt = PeriodicTable.GetPeriodicTable()
# Get a dictionary of atom counts {atomic_number: count}
atom_counts = rdMolDescriptors.GetAtomCounts(mol)

# Build the equation string for exact mass calculation
# to show "each number in the final equation"
equation_parts = []
print("Atom Counts in the final product:")
for atomic_num, count in sorted(atom_counts.items()):
    element_symbol = pt.GetElementSymbol(atomic_num)
    # Get the mass of the most common isotope for the element
    isotope_mass = pt.GetMostCommonIsotopeMass(atomic_num)
    
    # Print each number: atom count and isotope mass
    print(f"  - Element: {element_symbol}, Count: {count}, Isotope Mass: {isotope_mass:.4f}")
    equation_parts.append(f"({count} * {isotope_mass:.4f})")

# Calculate the final exact molecular weight
exact_mw = Descriptors.ExactMolWt(mol)

print("\nEquation for calculating exact mass:")
print(f"{' + '.join(equation_parts)} = {exact_mw:.4f} u")

print(f"\nFinal Exact Molecular Weight of Compound A: {exact_mw:.4f} u")