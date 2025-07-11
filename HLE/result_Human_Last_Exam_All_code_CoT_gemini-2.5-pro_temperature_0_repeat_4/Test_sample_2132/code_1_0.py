# Step 1: Define the parameters for cyclopentanecarboxylic acid.

# N: Total number of non-hydrogen atoms (6 carbons + 2 oxygens)
N = 8

# N_b: Number of non-hydrogen atoms in the cyclic part (5 carbons in the cyclopentane ring)
N_b = 5

# N_a: Number of non-hydrogen atoms in the acyclic part (1 carbon and 2 oxygens in the -COOH group)
N_a = 3

# Step 2: Calculate the Böttcher Molecular Complexity using the formula C = N * (N_b / N_a)^2.
bottcher_complexity = N * (N_b / N_a)**2

# Step 3: Print the final equation and the result.
print("The Böttcher Molecular Complexity (C) is calculated as follows:")
print(f"C = N * (N_b / N_a)^2")
print(f"C = {N} * ({N_b} / {N_a})^2")
print(f"C = {bottcher_complexity:.3f}")