# Define the parameters for the calculation
# d_c is the upper critical dimension for the n-vector model
d_c = 4
# d is the spatial dimension of interest
d = 3
# n is the number of components of the order parameter.
# We assume the Ising model universality class, where n=1.
n = 1

# Step 1: Calculate epsilon (ϵ)
epsilon = d_c - d

# Step 2: Calculate the value of the exponent alpha (α) using the first-order epsilon expansion formula
# Formula: alpha = (4 - n) / (2 * (n + 8)) * epsilon
numerator = 4 - n
denominator = 2 * (n + 8)
alpha = numerator / denominator * epsilon

# Step 3: Print the calculation step-by-step
print(f"The calculation for the specific heat exponent alpha (α) using the ϵ-expansion:")
print(f"Upper critical dimension, d_c = {d_c}")
print(f"Spatial dimension, d = {d}")
print(f"Number of order parameter components, n = {n} (Ising Model)")
print(f"Epsilon, ϵ = d_c - d = {d_c} - {d} = {epsilon}")
print(f"\nFormula: α = (4 - n) / (2 * (n + 8)) * ϵ")
print(f"Substituting the values:")
print(f"α = ({4} - {n}) / (2 * ({n} + 8)) * {epsilon}")
print(f"α = {numerator} / {denominator} * {epsilon}")
print(f"α = {alpha}")
print(f"\nTherefore, the scaling exponent α is {alpha:.4f} to first order in the epsilon expansion.")