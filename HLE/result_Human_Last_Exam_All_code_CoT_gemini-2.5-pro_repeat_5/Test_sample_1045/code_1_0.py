# Toluene molecule is C7H8
num_C = 7
num_H = 8

# Number of contracted Gaussian functions for each atom using the 6-31G basis set
# For Carbon (heavy atom): 1 core (1s) + 4 valence (inner 2s,2p) + 4 valence (outer 2s,2p) = 9
# Or more simply: 1(1s) + 2(2s) + 6(2px,2py,2pz) = 9
funcs_per_C = 9

# For Hydrogen: Valence shell is split into two parts (inner and outer) = 2
funcs_per_H = 2

# Calculate the total number of functions for the molecule
total_funcs_C = num_C * funcs_per_C
total_funcs_H = num_H * funcs_per_H
total_functions = total_funcs_C + total_funcs_H

# Print the breakdown of the calculation and the final answer
print(f"Molecule: Toluene (C{num_C}H{num_H})")
print("Basis Set: 6-31G")
print("-" * 30)
print(f"Functions per Carbon atom: {funcs_per_C}")
print(f"Functions per Hydrogen atom: {funcs_per_H}")
print("-" * 30)
print("Total number of contracted Gaussian functions is calculated as:")
print(f"({num_C} C atoms * {funcs_per_C} functions/C) + ({num_H} H atoms * {funcs_per_H} functions/H) = {total_functions}")
print(f"({total_funcs_C}) + ({total_funcs_H}) = {total_functions}")
print(f"\nFinal Answer: {total_functions}")
