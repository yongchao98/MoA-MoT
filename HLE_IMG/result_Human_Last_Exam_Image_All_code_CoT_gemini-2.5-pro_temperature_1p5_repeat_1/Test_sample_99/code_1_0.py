#
# This script calculates the molecular formula of the final product.
#

# --- Atom counts for the final product: N-benzyl-5-methyl-imidazo[2,1-b]thiazole-6-carboxamide ---

# Carbon count
# 5 (fused ring) + 1 (methyl) + 1 (carbonyl) + 7 (benzyl group)
num_C = 5 + 1 + 1 + 7

# Hydrogen count
# 2 (fused ring) + 3 (methyl) + 1 (amide NH) + 7 (benzyl group)
num_H = 2 + 3 + 1 + 7

# Nitrogen count
# 2 (fused ring) + 1 (amide)
num_N = 2 + 1

# Oxygen count
# 1 (amide carbonyl)
num_O = 1

# Sulfur count
# 1 (fused ring)
num_S = 1

# --- Construct and print the molecular formula ---

print("The molecular formula is composed of:")
print(f"Carbon (C): {num_C}")
print(f"Hydrogen (H): {num_H}")
print(f"Nitrogen (N): {num_N}")
print(f"Oxygen (O): {num_O}")
print(f"Sulfur (S): {num_S}")

molecular_formula = f"C{num_C}H{num_H}N{num_N}O{num_O}S{num_S}"
# O and S are often written without a subscript '1', let's adjust for that convention.
molecular_formula = f"C{num_C}H{num_H}N{num_N}OS"

print("\nThe final molecular formula is:")
print(molecular_formula)
