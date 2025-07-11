# This script calculates and displays the molecular formula of the final product.

# Step 1: Define the number of atoms for each element based on chemical analysis.
# The final product is N-benzyl-2-methylimidazo[2,1-b]thiazole-3-carboxamide.
# Let's count the atoms:
# Carbon (C): 5 (core) + 1 (methyl) + 1 (carbonyl) + 7 (benzyl) = 14
# Hydrogen (H): 2 (core) + 3 (methyl) + 1 (amide) + 7 (benzyl) = 13
# Nitrogen (N): 2 (core) + 1 (amide) = 3
# Oxygen (O): 1 (carbonyl)
# Sulfur (S): 1 (core)

num_C = 14
num_H = 13
num_N = 3
num_O = 1
num_S = 1

# Step 2: Construct the molecular formula string.
# The standard order is C, H, then other elements alphabetically.
# If an element count is 1, the number is typically omitted.
formula_parts = []
if num_C > 0:
    formula_parts.append(f"C{num_C}")
if num_H > 0:
    formula_parts.append(f"H{num_H}")
if num_N > 0:
    formula_parts.append(f"N{num_N}")
if num_O > 0:
    # Omitting the number '1' for single atoms as per convention
    formula_parts.append(f"O")
if num_S > 0:
    # Omitting the number '1' for single atoms as per convention
    formula_parts.append(f"S")

molecular_formula = "".join(formula_parts)

# Step 3: Print the final result.
# The prompt asks to output each number in the final equation/formula.
# We will print the formula with explicit numbers for clarity as requested.
print("The molecular formula of the product is composed of:")
print(f"Carbon (C): {num_C}")
print(f"Hydrogen (H): {num_H}")
print(f"Nitrogen (N): {num_N}")
print(f"Oxygen (O): {num_O}")
print(f"Sulfur (S): {num_S}")
print("\nFinal Molecular Formula:")
print(f"C{num_C}H{num_H}N{num_N}OS")
