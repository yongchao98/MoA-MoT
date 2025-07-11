# Define the atomic composition of reactants and products based on the reaction analysis.

# 1. Reactant cation: 1,3,6,8-tetramethoxy-9-phenylxanthenylium
reactant_cation = {'C': 23, 'H': 21, 'O': 5, 'N': 0}

# 2. Reagent for B: methyl-3-aminopropionate (NH2-CH2-CH2-COOCH3)
reagent_b = {'C': 4, 'H': 9, 'O': 2, 'N': 1}

# 3. Byproduct: water (H2O)
byproduct = {'C': 0, 'H': 2, 'O': 1, 'N': 0}

# 4. Counter-ion: tetrafluoroborate (BF4-)
counter_ion = {'B': 1, 'F': 4}

# --- Calculation ---
# The reaction is: ReactantCation + ReagentB -> ProductBCation + Water
# ProductBCation = ReactantCation + ReagentB - Water

print("Calculating the atomic composition of the cation of Compound B:")

c_cat = reactant_cation['C'] + reagent_b['C'] - byproduct['C']
print(f"Number of Carbon (C) atoms = {reactant_cation['C']} + {reagent_b['C']} - {byproduct['C']} = {c_cat}")

h_cat = reactant_cation['H'] + reagent_b['H'] - byproduct['H']
print(f"Number of Hydrogen (H) atoms = {reactant_cation['H']} + {reagent_b['H']} - {byproduct['H']} = {h_cat}")

n_cat = reactant_cation['N'] + reagent_b['N'] - byproduct['N']
print(f"Number of Nitrogen (N) atoms = {reactant_cation['N']} + {reagent_b['N']} - {byproduct['N']} = {n_cat}")

o_cat = reactant_cation['O'] + reagent_b['O'] - byproduct['O']
print(f"Number of Oxygen (O) atoms = {reactant_cation['O']} + {reagent_b['O']} - {byproduct['O']} = {o_cat}")

product_b_cation = {'C': c_cat, 'H': h_cat, 'N': n_cat, 'O': o_cat}

print(f"\nThe formula for the cation of Compound B is C{product_b_cation['C']}H{product_b_cation['H']}N{product_b_cation['N']}O{product_b_cation['O']}.")

# Compound B is the salt including the counter-ion BF4-.
product_b_salt = {
    'C': product_b_cation.get('C', 0),
    'H': product_b_cation.get('H', 0),
    'B': counter_ion.get('B', 0),
    'F': counter_ion.get('F', 0),
    'N': product_b_cation.get('N', 0),
    'O': product_b_cation.get('O', 0)
}

# Construct the final molecular formula string in standard order (C, H, then alphabetical)
# Omitting the number 1 for single atoms
formula_str = (
    f"C{product_b_salt['C']}"
    f"H{product_b_salt['H']}"
    f"B" + (str(product_b_salt['B']) if product_b_salt['B'] > 1 else "") +
    f"F{product_b_salt['F']}"
    f"N" + (str(product_b_salt['N']) if product_b_salt['N'] > 1 else "") +
    f"O{product_b_salt['O']}"
)

print("\nAdding the counter-ion BF4- gives the final molecular formula of Compound B.")
print(f"Final formula: {formula_str}")