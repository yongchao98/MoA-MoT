# Step 1: Define the number of backbone atoms for each monomer.
# Alanine is an alpha-amino acid with a backbone of (-NH-CH-CO-).
atoms_in_alanine = 3
# An epsilon-amino acid has a backbone of (-NH-(CH2)5-CO-).
atoms_in_epsilon_aa = 7

# Step 2: Calculate 'm', the number of backbone atoms in the repeating unit.
# The repeating unit consists of one alanine and one epsilon-amino acid.
m = atoms_in_alanine + atoms_in_epsilon_aa

# Step 3: Calculate 'n', the number of atoms in the most likely hydrogen-bonded ring.
# We analyze the i -> i+3 hydrogen bonding pattern, which is common and stable.
# The covalent loop for an i -> i+3 bond involves the backbones of residues i+1 and i+2.
# In an alternating sequence, these will be one alanine and one epsilon-amino acid.
# The ring size 'n' is the sum of atoms in the intervening backbones plus the two end atoms (C' and N).
# The order of intervening monomers doesn't change the sum.
n = 1 + atoms_in_alanine + atoms_in_epsilon_aa + 1

# Step 4: Output the results.
print("--- Calculation of Helical Pattern (m/n) ---")
print(f"Monomer 1 (Alanine) has {atoms_in_alanine} backbone atoms.")
print(f"Monomer 2 (epsilon-Amino Acid) has {atoms_in_epsilon_aa} backbone atoms.")
print(f"The number of atoms 'm' in the (Ala + εAA) repeating unit is: {atoms_in_alanine} + {atoms_in_epsilon_aa} = {m}")
print(f"The number of atoms 'n' in the most stable hydrogen-bonded ring (i -> i+3 pattern) is: 1 (C') + {atoms_in_alanine} (Ala) + {atoms_in_epsilon_aa} (εAA) + 1 (N) = {n}")
print(f"\nThe most likely helical pattern is m/n, which is {m}/{n}.")
