#
# Step 1: Explain the logic and define the products
#
# The invalid SMILES string likely represents a ketal that hydrolyzes in acid.
# The hydrolysis reaction is: Ketal + H2O -> Ketone + Diol.
# Based on the SMILES fragments, the most plausible products are Acetophenone and Isosorbide.
#
# Product 1: Acetophenone
# SMILES: CC(=O)c1ccccc1
#
# Product 2: Isosorbide
# SMILES: OC1C2OCC(O2)C1O
#
# We will now calculate their molar masses to find the heavier product.
#

try:
    from rdkit import Chem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
    import re
except ImportError:
    print("Please install RDKit to run this code: pip install rdkit")
    exit()

def calculate_and_print_molar_mass(name, smiles):
    """Calculates molar mass from a SMILES string and prints the detailed equation."""
    
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'O': 15.999,
        'N': 14.007,
        # Add other elements if needed
    }

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Error: Could not parse SMILES for {name}: {smiles}")
        return 0

    formula = CalcMolFormula(mol)
    
    # Parse the molecular formula string (e.g., "C8H8O")
    atom_counts = {element: 0 for element in atomic_weights}
    for element, count in re.findall(r'([A-Z][a-z]*)(\d*)', formula):
        if element in atom_counts:
            atom_counts[element] = int(count) if count else 1

    # Build the equation string and calculate the total mass
    total_mass = 0
    equation_parts = []
    for element, count in atom_counts.items():
        if count > 0:
            weight = atomic_weights[element]
            total_mass += count * weight
            equation_parts.append(f"{count} * {weight}")

    equation_str = " + ".join(equation_parts)
    print(f"For {name} ({formula}):")
    print(f"  Equation: Molar Mass = {equation_str}")
    print(f"  Result: {total_mass:.3f} g/mol\n")
    
    return total_mass, formula

# Define the two potential products
products = {
    "Acetophenone": "CC(=O)c1ccccc1",
    "Isosorbide": "OC1C2OCC(O2)C1O"
}

# Store results
results = {}

print("Calculating molar masses of the potential hydrolysis products...")
print("-" * 60)

for name, smiles in products.items():
    mass, formula = calculate_and_print_molar_mass(name, smiles)
    results[name] = {'smiles': smiles, 'mass': mass, 'formula': formula}

# Compare the molar masses to find the heavier product
product1_name = list(results.keys())[0]
product2_name = list(results.keys())[1]

product1 = results[product1_name]
product2 = results[product2_name]

if product1['mass'] > product2['mass']:
    heavier_product_name = product1_name
    heavier_product_smiles = product1['smiles']
else:
    heavier_product_name = product2_name
    heavier_product_smiles = product2['smiles']

print("Comparing the products...")
print(f"The product with the higher molar mass is {heavier_product_name}.")
print(f"Its SMILES string is: {heavier_product_smiles}")
print("-" * 60)
print("\nFinal Answer:")
print(heavier_product_smiles)

<<<OC1C2OCC(O2)C1O>>>