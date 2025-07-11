# This script calculates the molecular formula of the final product based on its structure.
# The reaction sequence leads to cis-4-amino-5-(trifluoromethyl)cyclohexane-1-carboxylic acid.
# We will count the atoms in this final structure.

# 1. Carbon (C) atom count
# There are 6 carbon atoms in the cyclohexane ring.
carbons_in_ring = 6
# There is 1 carbon atom in the carboxylic acid group (-COOH).
carbons_in_cooh = 1
# There is 1 carbon atom in the trifluoromethyl group (-CF3).
carbons_in_cf3 = 1
total_carbons = carbons_in_ring + carbons_in_cooh + carbons_in_cf3

# 2. Hydrogen (H) atom count
# On the cyclohexane ring, there are 9 hydrogen atoms.
# (CH-COOH, CH2, CH2, CH-NH2, CH-CF3, CH2) -> 1+2+2+1+1+2 = 9
hydrogens_on_ring = 9
# There is 1 hydrogen in the carboxylic acid group (-COOH).
hydrogens_in_cooh = 1
# There are 2 hydrogens in the amino group (-NH2).
hydrogens_in_nh2 = 2
total_hydrogens = hydrogens_on_ring + hydrogens_in_cooh + hydrogens_in_nh2

# 3. Fluorine (F) atom count
# There are 3 fluorine atoms in the trifluoromethyl group (-CF3).
total_fluorines = 3

# 4. Nitrogen (N) atom count
# There is 1 nitrogen atom in the amino group (-NH2).
total_nitrogens = 1

# 5. Oxygen (O) atom count
# There are 2 oxygen atoms in the carboxylic acid group (-COOH).
total_oxygens = 2

# As requested, printing the calculation for each element in an equation-like format.
print("Calculation of the molecular formula based on the final product's structure:")
print(f"Number of Carbon atoms (C) = {carbons_in_ring} (ring) + {carbons_in_cooh} (COOH) + {carbons_in_cf3} (CF3) = {total_carbons}")
print(f"Number of Hydrogen atoms (H) = {hydrogens_on_ring} (ring) + {hydrogens_in_cooh} (COOH) + {hydrogens_in_nh2} (NH2) = {total_hydrogens}")
print(f"Number of Fluorine atoms (F) = {total_fluorines} (CF3) = {total_fluorines}")
print(f"Number of Nitrogen atoms (N) = {total_nitrogens} (NH2) = {total_nitrogens}")
print(f"Number of Oxygen atoms (O) = {total_oxygens} (COOH) = {total_oxygens}")

# Construct and print the final molecular formula in the standard order (C, H, then alphabetical).
molecular_formula = f"C{total_carbons}H{total_hydrogens}F{total_fluorines}N{total_nitrogens}O{total_oxygens}"
print(f"\nThe final molecular formula is: {molecular_formula}")