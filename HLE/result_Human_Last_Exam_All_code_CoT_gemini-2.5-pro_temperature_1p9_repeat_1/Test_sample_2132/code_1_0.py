import math

# This script calculates the Böttcher Molecular Complexity for the product of the
# Favorskii rearrangement of 2-chlorocyclohexanone.

# Step 1: Identify the product and its properties
product_name = "cyclopentanecarboxylic acid"
# N_A: Number of non-hydrogen atoms (6 Carbons + 2 Oxygens)
N_A = 8
# N_B: Number of bonds between non-hydrogen atoms (5 in ring + 1 link + 2 in COOH)
N_B = 8
# N_stereo: Number of stereoisomers (1 chiral center -> 2^1 stereoisomers)
N_stereo = 2

# Print an explanation of the process
print(f"The product is {product_name}.")
print("To find its Böttcher Molecular Complexity, we use the formula:")
print("Complexity = (N_A * N_B * log₂(N_A * N_B)) / N_stereo\n")

print("Calculating the parameters:")
print(f"Number of non-hydrogen atoms (N_A) = {N_A}")
print(f"Number of bonds between non-hydrogen atoms (N_B) = {N_B}")
print(f"Number of stereoisomers (N_stereo) = {N_stereo}\n")

# Step 2: Perform the calculation
N_A_times_N_B = N_A * N_B
log_value = math.log2(N_A_times_N_B)
complexity = (N_A_times_N_B * log_value) / N_stereo

# Step 3: Display the final equation and result
# We use int() to show clean integer values where appropriate
print("Final equation with values:")
print(f"({N_A} * {N_B} * log₂({N_A_times_N_B})) / {N_stereo} = ({N_A_times_N_B} * {int(log_value)}) / {N_stereo} = {int(complexity)}")
