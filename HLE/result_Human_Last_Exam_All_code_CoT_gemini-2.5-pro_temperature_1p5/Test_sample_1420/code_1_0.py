#
# Plan:
# 1. Define the parameters for the epsilon expansion calculation.
#    - N: The number of components of the order parameter (N=1 for the Ising model).
#    - d: The spatial dimension of the system.
#    - d_c: The upper critical dimension for this model.
# 2. Calculate the expansion parameter, epsilon = d_c - d.
# 3. Use the first-order expansion formula for the specific heat exponent alpha.
#    alpha = (4 - N) / (2 * (N + 8)) * epsilon
# 4. Substitute the parameters into the formula and print the step-by-step calculation.
# 5. Print the final result.
#

# Step 1: Define parameters
N = 1   # For the Ising universality class
d = 3   # Spatial dimension given in the problem
d_c = 4 # Upper critical dimension for the O(N) model

# Step 2: Calculate epsilon
epsilon = d_c - d

# Step 3 & 4: Use the formula and calculate alpha
# The formula for alpha to first order in epsilon is:
# alpha = (4 - N) / (2 * (N + 8)) * epsilon

numerator = 4 - N
denominator_part1 = 2
denominator_part2 = N + 8
denominator = denominator_part1 * denominator_part2

alpha = (numerator / denominator) * epsilon

# Step 5: Print the results and the calculation steps
print(f"Calculating the specific heat exponent α using the ε-expansion for d={d}.")
print(f"The upper critical dimension is d_c = {d_c}.")
print(f"The expansion parameter ε = d_c - d = {d_c} - {d} = {epsilon}.")
print("Using the O(N) model with N=1 (Ising class).")
print("\nThe first-order formula for α is: α = (4 - N) / (2 * (N + 8)) * ε")
print("\nSubstituting the values:")
print(f"α = ({4} - {N}) / ({denominator_part1} * ({N} + {8})) * {epsilon}")
print(f"α = {numerator} / ({denominator_part1} * {denominator_part2}) * {epsilon}")
print(f"α = {numerator} / {denominator} * {epsilon}")
print(f"α = {alpha}")
