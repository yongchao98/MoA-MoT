import math

# Step 1: Define the product of the reaction
product_name = "Cyclopentanecarboxylic acid"
molecular_formula = "C6H10O2"

print(f"The product of the Favorskii rearrangement of 2-chlorocyclohexanone is {product_name}.")
print("-" * 30)

# Step 2: Count the parameters for the Böttcher Molecular Complexity Index.
# For Cyclopentanecarboxylic acid:
# The structure is a 5-membered carbon ring attached to a -COOH group.

# N_a: Number of non-hydrogen atoms (6 carbons + 2 oxygens)
N_a = 6 + 2

# N_b: Number of bonds between non-hydrogen atoms.
# 5 bonds in the cyclopentane ring (C-C)
# 1 bond connecting the ring to the carboxyl carbon (C-C)
# 2 bonds within the carboxyl group (C=O and C-O)
# We count connections, so a double bond is 1 connection.
N_b = 5 + 1 + 2

# N_c: Cyclomatic number (number of rings). For a single connected molecule, N_c = N_b - N_a + 1.
N_c = N_b - N_a + 1

print("Calculating the parameters for the Böttcher Molecular Complexity Index:")
print(f"Number of non-hydrogen atoms (N_a) = {N_a}")
print(f"Number of bonds between non-hydrogen atoms (N_b) = {N_b}")
print(f"Cyclomatic number (number of rings, N_c) = {N_c}")
print("-" * 30)

# Step 3: Calculate the Böttcher Molecular Complexity Index (MCI)
# Formula: MCI = (N_a * N_b) / (N_a + N_b - N_c)
numerator = N_a * N_b
denominator = N_a + N_b - N_c
mci_result = numerator / denominator

print("The formula for the index is: MCI = (N_a * N_b) / (N_a + N_b - N_c)")
print(f"Plugging in the values: MCI = ({N_a} * {N_b}) / ({N_a} + {N_b} - {N_c})")
print(f"MCI = {numerator} / {denominator}")
print(f"The final Böttcher Molecular Complexity Index is: {mci_result}")
