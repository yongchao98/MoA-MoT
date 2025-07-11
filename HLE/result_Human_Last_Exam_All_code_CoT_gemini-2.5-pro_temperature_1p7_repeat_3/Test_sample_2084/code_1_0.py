# Import necessary libraries from RDKit
from rdkit import Chem
from rdkit.Chem import Descriptors
import re

def calculate_and_print_mass(name, smiles):
    """
    Calculates the molar mass of a molecule from its SMILES string,
    printing the detailed breakdown of the calculation.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Error: Could not parse SMILES for {name}")
        return 0, None

    # Define atomic weights
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'O': 15.999,
        'N': 14.007,
        'F': 18.998,
        'Cl': 35.453,
        # Add other elements if needed
    }

    # Get the molecular formula and parse it
    formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol)
    
    # Use regex to find all element-count pairs (e.g., C2, H4, O)
    atom_counts = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
    
    total_mass = 0
    print(f"\nCalculating molar mass for {name} ({formula}):")

    for element, count_str in atom_counts:
        count = int(count_str) if count_str else 1
        weight = atomic_weights.get(element, 0)
        mass_contrib = count * weight
        total_mass += mass_contrib
        # Output the calculation for this element
        print(f"  {element}: {count} * {weight:.3f} = {mass_contrib:.3f}")

    print(f"Total Molar Mass for {name} = {total_mass:.3f} g/mol")
    return total_mass, mol

# The two products from the hydrolysis reaction
product1_name = "Acetone"
product1_smiles = "CC(=O)C"

product2_name = "1-phenyl-1,2-ethanediol"
product2_smiles = "OCC(O)c1ccccc1"

# Calculate masses and get mol objects
mass1, mol1 = calculate_and_print_mass(product1_name, product1_smiles)
mass2, mol2 = calculate_and_print_mass(product2_name, product2_smiles)

# Compare the masses and determine the product with the higher molar mass
print("\n--- Comparison ---")
if mass1 > mass2:
    print(f"{product1_name} has a higher molar mass.")
    # Canonical SMILES are generally preferred for output
    final_answer_smiles = Chem.MolToSmiles(mol1)
else:
    print(f"{product2_name} has a higher molar mass.")
    final_answer_smiles = Chem.MolToSmiles(mol2)

print("\nThe SMILES string for the product with the higher molar mass is:")
print(final_answer_smiles)
