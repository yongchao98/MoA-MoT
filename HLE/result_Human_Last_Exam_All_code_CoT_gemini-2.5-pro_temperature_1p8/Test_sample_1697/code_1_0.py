# The user needs to have rdkit installed. If not, they can install it via pip:
# pip install rdkit

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import re

print("--- Step-by-Step Chemical Analysis ---")
print("1. Starting Material: N,N-diethyl-3-dimethylaminobenzamide.")
print("2. Reagents: a) sec-BuLi/TMEDA, a powerful base for directed lithiation. b) Methyl iodide (CH3I), a methylating agent.")
print("3. Mechanism (Directed Ortho-Metalation):")
print("   - The amide (-CONEt2) and dimethylamino (-NMe2) groups are both excellent directing groups.")
print("   - The proton at the C2 position of the benzene ring is ortho to BOTH groups, making it the most acidic proton.")
print("   - sec-BuLi removes this C2 proton, creating a highly reactive aryllithium intermediate.")
print("4. Final Step (Electrophilic Quench):")
print("   - The aryllithium attacks the methyl group of methyl iodide, adding a methyl substituent to the C2 position.")
print("-" * 35)

# The product is N,N-diethyl-2-methyl-3-dimethylaminobenzamide.
# We can represent it with a SMILES string.
product_smiles = "O=C(N(CC)CC)c1c(C)c(N(C)C)ccc1"
product_mol = Chem.MolFromSmiles(product_smiles)

if product_mol:
    product_name = "N,N-diethyl-2-methyl-3-dimethylaminobenzamide"
    molecular_formula = rdMolDescriptors.CalcMolFormula(product_mol)
    molecular_weight = Descriptors.MolWt(product_mol)

    print("--- Final Product Information ---")
    print(f"Compound Name: {product_name}")
    print(f"SMILES String: {product_smiles}")
    print(f"Molecular Formula: {molecular_formula}")
    print(f"Molecular Weight: {molecular_weight:.2f}")
    print("-" * 35)

    # The following section interprets "output each number in the final equation"
    # by showing the count of each atom in the final product's formula.
    print("--- Product's Atomic Composition (The Final Equation) ---")
    print(f"The final product's molecular formula is {molecular_formula}.")
    
    # Use regular expression to find element symbols and their counts from the formula string
    atomic_composition = re.findall(r'([A-Z][a-z]*)(\d*)', molecular_formula)
    
    equation_parts = []
    for element, count in atomic_composition:
        # If the count is not specified in the formula, it defaults to 1
        number = int(count) if count else 1
        equation_parts.append(f"Number of {element} atoms: {number}")

    # Print each part of the "equation"
    for part in equation_parts:
        print(part)

else:
    print("Error: Could not process the chemical structure from the SMILES string.")

<<<N,N-diethyl-2-methyl-3-dimethylaminobenzamide>>>