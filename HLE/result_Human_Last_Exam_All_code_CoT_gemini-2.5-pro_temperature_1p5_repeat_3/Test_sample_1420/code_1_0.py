# Parameters for the calculation
# The spatial dimension of the system
d = 3
# The upper critical dimension for the n-vector model
d_c = 4
# The number of components of the order parameter. We assume the Ising model (n=1).
n = 1

# 1. Calculate epsilon (ϵ)
epsilon = d_c - d

# 2. Calculate the specific heat exponent alpha (α) to first order in epsilon
# The formula for alpha is (4 - n) / (2 * (n + 8)) * epsilon
alpha_numerator = 4 - n
alpha_denominator = 2 * (n + 8)
alpha = (alpha_numerator / alpha_denominator) * epsilon

# 3. Print the results step-by-step
print(f"The calculation is for the {n}-vector model (Ising universality class).")
print(f"Spatial dimension d = {d}")
print(f"Upper critical dimension d_c = {d_c}")
print(f"Epsilon (ϵ = d_c - d) = {d_c} - {d} = {epsilon}")
print("\nThe formula for the specific heat exponent α is:")
print("α = (4 - n) / (2 * (n + 8)) * ϵ")
print("\nSubstituting the values:")
print(f"α = ({4} - {n}) / (2 * ({n} + {8})) * {epsilon}")
print(f"α = {alpha_numerator} / {alpha_denominator} * {epsilon}")
print(f"α = {alpha_numerator / alpha_denominator} * {epsilon}")

# Final result
print("\nThe scaling exponent α for the specific heat is:")
print(alpha)