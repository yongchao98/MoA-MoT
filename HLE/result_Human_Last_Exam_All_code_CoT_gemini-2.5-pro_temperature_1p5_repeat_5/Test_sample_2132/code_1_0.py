import math

# Step 1: Define the properties of the product molecule, cyclopentanecarboxylic acid.

# N is the number of non-hydrogen atoms.
# Cyclopentanecarboxylic acid (C6H10O2) has 6 Carbon atoms and 2 Oxygen atoms.
N = 6 + 2

# E is the number of bonds between non-hydrogen atoms.
# - 5 C-C bonds in the cyclopentane ring.
# - 1 C-C bond from the ring to the carboxyl group.
# - 1 C-O single bond in the carboxyl group.
# - 1 C=O double bond in the carboxyl group (counted as one connection).
E = 5 + 1 + 1 + 1

# Step 2: Apply the Böttcher Molecular Complexity formula.
# M = N * (E / (N - 1))
if N > 1:
    bottcher_complexity = N * (E / (N - 1))
else:
    # The formula is undefined for N <= 1.
    bottcher_complexity = 0

# Step 3: Print the explanation and the final result.
print("The product of the Favorskii rearrangement of 2-chlorocyclohexanone is cyclopentanecarboxylic acid.")
print("The Böttcher Molecular Complexity (M) is calculated using the formula: M = N * (E / (N - 1))")
print(f"Number of non-hydrogen atoms (N) = {N}")
print(f"Number of bonds between non-hydrogen atoms (E) = {E}")
print("\nCalculating the complexity:")
# The problem asks to output the final equation with the numbers.
print(f"M = {N} * ({E} / ({N} - 1))")
print(f"M = {N} * ({E} / {N-1})")
print(f"M = {N*E} / {N-1}")
print(f"M = {bottcher_complexity}")
