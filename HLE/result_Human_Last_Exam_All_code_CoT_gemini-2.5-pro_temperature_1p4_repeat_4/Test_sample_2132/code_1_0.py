#
# This script calculates the Böttcher Molecular Complexity for the product of the
# Favorskii rearrangement of 2-chlorocyclohexanone.
#

# Step 1: Identify the reaction product.
# The Favorskii rearrangement of 2-chlorocyclohexanone produces cyclopentanecarboxylic acid.
product_name = "cyclopentanecarboxylic acid"
product_formula = "C6H10O2"

print(f"The product of the Favorskii rearrangement of 2-chlorocyclohexanone is {product_name}.")
print(f"The chemical formula of the product is {product_formula}.")
print("-" * 20)

# Step 2: Define the Böttcher Molecular Complexity formula.
# The formula is BC = Na * Nb / Np
# Na = total number of atoms
# Nb = total number of covalent bonds (double bonds count as 2)
# Np = number of connected components (1 for a single molecule)

print("The Böttcher Molecular Complexity is calculated as: Na * Nb / Np")
print("-" * 20)

# Step 3: Calculate Na, Nb, and Np for the product (C6H10O2).

# Na: Number of atoms
num_C = 6
num_H = 10
num_O = 2
Na = num_C + num_H + num_O
print(f"Calculating the number of atoms (Na):")
print(f"  - Carbon atoms: {num_C}")
print(f"  - Hydrogen atoms: {num_H}")
print(f"  - Oxygen atoms: {num_O}")
print(f"Total Na = {num_C} + {num_H} + {num_O} = {Na}")
print("-" * 20)


# Nb: Number of covalent bonds
# Calculated as (sum of valences of all atoms) / 2
# Valence: C=4, H=1, O=2
sum_of_valences = (num_C * 4) + (num_H * 1) + (num_O * 2)
Nb = sum_of_valences // 2
print(f"Calculating the number of bonds (Nb):")
print(f"  - Sum of valences = ({num_C} * 4) + ({num_H} * 1) + ({num_O} * 2) = {sum_of_valences}")
print(f"Total Nb = {sum_of_valences} / 2 = {Nb}")
print("-" * 20)


# Np: Number of connected components
Np = 1
print(f"Calculating the number of components (Np):")
print(f"  - The product is a single molecule, so Np = {Np}")
print("-" * 20)

# Step 4: Calculate the final Böttcher Molecular Complexity.
complexity = Na * Nb / Np

print("Final Calculation:")
print(f"The final equation is: {Na} * {Nb} / {Np} = {int(complexity)}")