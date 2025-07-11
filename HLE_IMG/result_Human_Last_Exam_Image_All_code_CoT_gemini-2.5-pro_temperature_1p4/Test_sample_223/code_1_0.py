# Define the atomic composition of the reactants and byproducts
start_cation = {'C': 25, 'H': 27, 'N': 0, 'O': 6}
amine_B = {'C': 4, 'H': 9, 'N': 1, 'O': 2}
byproducts = {'C': 2, 'H': 8, 'N': 0, 'O': 2}

# Calculate the atomic composition of the product cation B
product_B_cation = {}
product_B_cation['C'] = start_cation['C'] + amine_B['C'] - byproducts['C']
product_B_cation['H'] = start_cation['H'] + amine_B['H'] - byproducts['H']
product_B_cation['N'] = start_cation['N'] + amine_B['N'] - byproducts['N']
product_B_cation['O'] = start_cation['O'] + amine_B['O'] - byproducts['O']

# Print the calculation steps
print("Calculating the molecular formula for the cation of compound B:")
print(f"Carbon atoms: {start_cation['C']} (start) + {amine_B['C']} (amine) - {byproducts['C']} (byproduct) = {product_B_cation['C']}")
print(f"Hydrogen atoms: {start_cation['H']} (start) + {amine_B['H']} (amine) - {byproducts['H']} (byproduct) = {product_B_cation['H']}")
print(f"Nitrogen atoms: {start_cation['N']} (start) + {amine_B['N']} (amine) - {byproducts['N']} (byproduct) = {product_B_cation['N']}")
print(f"Oxygen atoms: {start_cation['O']} (start) + {amine_B['O']} (amine) - {byproducts['O']} (byproduct) = {product_B_cation['O']}")

# Format and print the final molecular formula
formula = f"C{product_B_cation['C']}H{product_B_cation['H']}N{product_B_cation['N'] if product_B_cation['N'] > 1 else ''}O{product_B_cation['O']}"
print(f"\nThe molecular formula of the cation of compound B is: {formula}")
