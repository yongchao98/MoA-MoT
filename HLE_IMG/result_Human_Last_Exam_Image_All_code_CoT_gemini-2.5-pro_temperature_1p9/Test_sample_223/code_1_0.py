# Define the atomic composition of the reactants and byproducts.

# 1. Starting material cation (1,8-dimethoxy-9-(2,6-dimethoxyphenyl)xanthenylium)
sm_cation = {'C': 23, 'H': 21, 'N': 0, 'O': 5}

# 2. Reactant amine (methyl-3-aminopropionate)
amine = {'C': 4, 'H': 9, 'N': 1, 'O': 2}

# 3. Eliminated byproduct (water)
water = {'C': 0, 'H': 2, 'N': 0, 'O': 1}

# Calculate the atomic composition of the product cation B+
# Formula: Product = Reactant1 + Reactant2 - Byproduct
product_b_cation = {}

product_b_cation['C'] = sm_cation['C'] + amine['C'] - water['C']
product_b_cation['H'] = sm_cation['H'] + amine['H'] - water['H']
product_b_cation['N'] = sm_cation['N'] + amine['N'] - water['N']
product_b_cation['O'] = sm_cation['O'] + amine['O'] - water['O']

# Print the calculation steps for each element
print("Calculation for the molecular formula of compound B cation:")
print(f"Carbon (C): {sm_cation['C']} + {amine['C']} - {water['C']} = {product_b_cation['C']}")
print(f"Hydrogen (H): {sm_cation['H']} + {amine['H']} - {water['H']} = {product_b_cation['H']}")
print(f"Nitrogen (N): {sm_cation['N']} + {amine['N']} - {water['N']} = {product_b_cation['N']}")
print(f"Oxygen (O): {sm_cation['O']} + {amine['O']} - {water['O']} = {product_b_cation['O']}")
print("\n")

# Format and print the final molecular formula for the organic part of compound B
molecular_formula_B = f"C{product_b_cation['C']}H{product_b_cation['H']}N{product_b_cation['N']}O{product_b_cation['O']}"
print(f"The molecular formula of the organic cation in compound B is: {molecular_formula_B}")
