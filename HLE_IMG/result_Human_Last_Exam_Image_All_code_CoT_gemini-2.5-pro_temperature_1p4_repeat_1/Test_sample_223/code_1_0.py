import collections

def get_molecular_formula(atom_counts):
    """Formats a dictionary of atom counts into a molecular formula string."""
    formula_parts = []
    # Standard order: C, H, then alphabetical
    if 'C' in atom_counts:
        formula_parts.append(f"C{atom_counts['C']}")
    if 'H' in atom_counts:
        formula_parts.append(f"H{atom_counts['H']}")
    
    other_elements = sorted([el for el in atom_counts if el not in ['C', 'H']])
    for element in other_elements:
        count = atom_counts[element]
        # Do not show subscript if count is 1
        if count > 1:
            formula_parts.append(f"{element}{count}")
        else:
            formula_parts.append(element)
            
    return "".join(formula_parts)

# Step 1: Define the molecular formula of the starting cation as a dictionary of atom counts.
# The starting cation is 1,8-dimethoxy-9-(2,6-dimethoxyphenyl)xanthenylium, [C23H21O5]+.
starting_cation = collections.Counter({'C': 23, 'H': 21, 'O': 5})

# Step 2: Define the molecular formula of the reactant, methyl-3-aminopropionate.
# Structure: H2N-CH2-CH2-COOCH3, which is C4H9NO2.
reactant = collections.Counter({'C': 4, 'H': 9, 'N': 1, 'O': 2})

# Step 3: Define the byproduct of the reaction, which is water (H2O).
byproduct = collections.Counter({'H': 2, 'O': 1})

# Step 4: Calculate the molecular formula of the product cation (B)
# Product Cation = Starting Cation + Reactant - Byproduct
product_b_cation = starting_cation + reactant - byproduct

# Step 5: Print the explanation and the final molecular formula.
print("To find the molecular formula of the cation of compound B, we sum the atoms of the reactants and subtract the atoms of the byproduct.")
print("\nCalculation breakdown:")
print(f"Carbon (C): {starting_cation['C']} (from cation) + {reactant['C']} (from amine) = {product_b_cation['C']}")
print(f"Hydrogen (H): {starting_cation['H']} (from cation) + {reactant['H']} (from amine) - {byproduct['H']} (from H2O) = {product_b_cation['H']}")
print(f"Nitrogen (N): {starting_cation.get('N', 0)} (from cation) + {reactant['N']} (from amine) = {product_b_cation['N']}")
print(f"Oxygen (O): {starting_cation['O']} (from cation) + {reactant['O']} (from amine) - {byproduct['O']} (from H2O) = {product_b_cation['O']}")

# Format the final formula string
molecular_formula_b = get_molecular_formula(product_b_cation)

print(f"\nThe resulting molecular formula for the cation of compound B is {molecular_formula_b}.")
