# Step 1: Define the molecular formula of compound A based on its structure.
# Compound A is N-propyl-tetramethoxyacridinium cation.
# Core: C13 H5 N (Tetrasubstituted Acridinium C-skeleton)
# Substituents: 4x(-OCH3), 1x(-CH2CH2CH3)
# C = 13 (core) + 4*1 (methoxy) + 3 (propyl) = 20
# H = 5 (core) + 4*3 (methoxy) + 7 (propyl) = 24
# N = 1 (core)
# O = 4*1 (methoxy) = 4
compound_A = {'C': 20, 'H': 24, 'N': 1, 'O': 4}

# Step 2: Define the formula of the N-substituent group in A (n-propyl).
# R_A = -CH2CH2CH3
propyl_group = {'C': 3, 'H': 7, 'N': 0, 'O': 0}

# Step 3: Define the formula of the N-substituent group in B (from methyl-3-aminopropionate).
# Amine: H2N-CH2CH2COOCH3
# Group R_B = -CH2CH2COOCH3
# C = 1 + 1 + 1 + 1 = 4
# H = 2 + 2 + 3 = 7
# O = 2
new_group = {'C': 4, 'H': 7, 'N': 0, 'O': 2}

# Step 4: Calculate the molecular formula of compound B by swapping the groups.
# Formula of B = Formula of A - Formula of propyl_group + Formula of new_group
c_b = compound_A['C'] - propyl_group['C'] + new_group['C']
h_b = compound_A['H'] - propyl_group['H'] + new_group['H']
n_b = compound_A['N'] - propyl_group['N'] + new_group['N']
o_b = compound_A['O'] - propyl_group['O'] + new_group['O']

print("Calculating the molecular formula for compound B:")
print(f"Carbon atoms: {compound_A['C']} (from A) - {propyl_group['C']} (from propyl) + {new_group['C']} (from new group) = {c_b}")
print(f"Hydrogen atoms: {compound_A['H']} (from A) - {propyl_group['H']} (from propyl) + {new_group['H']} (from new group) = {h_b}")
print(f"Nitrogen atoms: {compound_A['N']} (from A) - {propyl_group['N']} (from propyl) + {new_group['N']} (from new group) = {n_b}")
print(f"Oxygen atoms: {compound_A['O']} (from A) - {propyl_group['O']} (from propyl) + {new_group['O']} (from new group) = {o_b}")

print(f"\nThe molecular formula of compound B is C{c_b}H{h_b}NO{o_b}.")
