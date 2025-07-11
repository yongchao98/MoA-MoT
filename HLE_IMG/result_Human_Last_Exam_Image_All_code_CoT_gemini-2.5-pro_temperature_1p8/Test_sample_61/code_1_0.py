# Molecular formula of Compound 1 (Methyl 1-oxo-2,3-dihydro-1H-indene-2-carboxylate)
compound_1 = {'C': 11, 'H': 10, 'O': 3}

# Atoms in the group being removed: Methoxycarbonyl (-COOMe)
group_removed = {'C': 2, 'H': 3, 'O': 2}

# Atoms in the group being added: Benzyl (-C7H7) from Compound 2
group_added = {'C': 7, 'H': 7, 'O': 0}

# Calculate the molecular formula of the final product A
product_A = {}
product_A['C'] = compound_1['C'] - group_removed['C'] + group_added['C']
product_A['H'] = compound_1['H'] - group_removed['H'] + group_added['H']
product_A['O'] = compound_1['O'] - group_removed['O'] + group_added['O']

# Print the step-by-step calculation
print("This reaction involves alkylation, followed by saponification and decarboxylation.")
print("The net transformation is the replacement of the methoxycarbonyl (-COOMe) group with a benzyl group.")
print("\nStarting with Compound 1 (C11H10O3).")
print(f"Group removed (-COOMe): C={group_removed['C']}, H={group_removed['H']}, O={group_removed['O']}")
print(f"Group added (benzyl): C={group_added['C']}, H={group_added['H']}, O={group_added['O']}")

print("\n--- Calculation for Product A ---")
print(f"Carbon atoms (C): {compound_1['C']} - {group_removed['C']} + {group_added['C']} = {product_A['C']}")
print(f"Hydrogen atoms (H): {compound_1['H']} - {group_removed['H']} + {group_added['H']} = {product_A['H']}")
print(f"Oxygen atoms (O): {compound_1['O']} - {group_removed['O']} + {group_added['O']} = {product_A['O']}")

# Format the final formula as a string
formula_A_str = f"C{product_A['C']}H{product_A['H']}O{product_A['O'] if product_A['O'] > 1 else ''}"

print(f"\nThe molecular formula of compound A is: {formula_A_str}")
