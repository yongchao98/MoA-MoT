# This script calculates the molecular formula of the product A.

# 1. Define the molecular formula of the starting material, compound 1.
c_start = 11
h_start = 10
o_start = 3

print("The calculation for the molecular formula of compound A is as follows:")

# 2. Calculate the final number of Carbon atoms.
# We add 7 carbons from the benzyl group and remove 2 from the -COOMe group.
c_benzyl = 7
c_coome = 2
c_final = c_start + c_benzyl - c_coome
print(f"\nCarbon atoms = {c_start} (from reactant 1) + {c_benzyl} (from benzyl group) - {c_coome} (from -COOMe group) = {c_final}")

# 3. Calculate the final number of Hydrogen atoms.
# We remove 1 acidic H, add 7 from the benzyl group, remove 3 from the -COOMe group, and add 1 back.
h_acidic = 1
h_benzyl = 7
h_coome = 3
h_replacement = 1
h_final = h_start - h_acidic + h_benzyl - h_coome + h_replacement
print(f"Hydrogen atoms = {h_start} (from reactant 1) - {h_acidic} (acidic H) + {h_benzyl} (from benzyl group) - {h_coome} (from -COOMe group) + {h_replacement} (replacement H) = {h_final}")

# 4. Calculate the final number of Oxygen atoms.
# We remove 2 oxygens from the -COOMe group.
o_coome = 2
o_final = o_start - o_coome
print(f"Oxygen atoms = {o_start} (from reactant 1) - {o_coome} (from -COOMe group) = {o_final}")

# 5. Assemble and print the final molecular formula.
# For elements with a count of 1, the number is omitted in the formula.
final_formula = f"C{c_final}H{h_final}O"
print(f"\nThe molecular formula of compound A is {final_formula}.")
