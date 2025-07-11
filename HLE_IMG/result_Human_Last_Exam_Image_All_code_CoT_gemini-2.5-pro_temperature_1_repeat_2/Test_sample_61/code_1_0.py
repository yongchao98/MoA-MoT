# Function to print the molecular formula from a dictionary of atom counts
def print_formula(atoms, name):
    formula_str = f"The molecular formula of {name} is "
    # Order: C, H, then alphabetically for other elements
    formula_str += f"C{atoms.get('C', 0)}"
    formula_str += f"H{atoms.get('H', 0)}"
    other_elements = sorted([el for el in atoms if el not in ['C', 'H']])
    for el in other_elements:
        if atoms[el] > 1:
            formula_str += f"{el}{atoms[el]}"
        elif atoms[el] == 1:
            formula_str += el
    print(formula_str)

# Step 1: Define the atoms in the starting material and the groups involved in the reaction.
compound_1 = {'C': 11, 'H': 10, 'O': 3}
benzyl_group = {'C': 7, 'H': 7}
coome_group = {'C': 2, 'H': 3, 'O': 2}

print("Step-by-step calculation for the molecular formula of Compound A:")
print("-" * 60)
print(f"1. Start with Compound 1: C{compound_1['C']}H{compound_1['H']}O{compound_1['O']}")

# Step 2: Calculate the formula of the alkylated intermediate.
# This involves removing one acidic hydrogen and adding a benzyl group.
alkylated_intermediate = {}
alkylated_intermediate['C'] = compound_1['C'] + benzyl_group['C']
alkylated_intermediate['H'] = compound_1['H'] - 1 + benzyl_group['H']
alkylated_intermediate['O'] = compound_1['O']

print(f"2. Alkylation: Replace one H with a benzyl group (C{benzyl_group['C']}H{benzyl_group['H']}).")
print(f"   Formula of intermediate: C{alkylated_intermediate['C']}H{alkylated_intermediate['H']}O{alkylated_intermediate['O']}")


# Step 3: Calculate the formula of the final product A.
# This involves saponification and decarboxylation, which replaces the -COOMe group with a hydrogen atom.
product_A = {}
product_A['C'] = alkylated_intermediate['C'] - coome_group['C']
product_A['H'] = alkylated_intermediate['H'] - coome_group['H'] + 1
product_A['O'] = alkylated_intermediate['O'] - coome_group['O']

print(f"3. Saponification & Decarboxylation: Replace -COOMe group (C{coome_group['C']}H{coome_group['H']}O{coome_group['O']}) with H.")

print("\nFinal calculation:")
print(f"Carbon atoms   = {compound_1['C']} (from Cpd 1) - 1 (H removed) + {benzyl_group['C']} (benzyl added) - {coome_group['C']} (COOMe removed) = {product_A['C']}") # This calculation is confusing. Let's show the final step-by-step calculation
print(f"Carbon atoms   = {compound_1['C']} + {benzyl_group['C']} - {coome_group['C']} = {product_A['C']}")
print(f"Hydrogen atoms = {compound_1['H']} - 1 (H removed) + {benzyl_group['H']} (benzyl added) - {coome_group['H']} (COOMe removed) + 1 (H added) = {product_A['H']}")
print(f"Oxygen atoms   = {compound_1['O']} - {coome_group['O']} (COOMe removed) = {product_A['O']}")
print("-" * 60)

# Print the final molecular formula of product A
print_formula(product_A, "compound A")