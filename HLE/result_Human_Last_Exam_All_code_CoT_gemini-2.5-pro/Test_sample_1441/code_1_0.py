# Define the parameters for the calculation
# d: spatial dimensionality
# N: number of components of the scalar field (N=1 for the Ising universality class)
d = 3
N = 1

# The upper critical dimension for this model
d_c = 4

# The epsilon expansion parameter
epsilon = d_c - d

# The formula for the specific heat exponent alpha to one-loop order is:
# alpha = (4 - N) / (N + 8) * epsilon
numerator = 4 - N
denominator = N + 8
alpha = (numerator / denominator) * epsilon

# Print the context and the step-by-step calculation
print(f"Calculating the specific heat exponent 'α' for a scalar field theory (N={N}) in d={d} dimensions.")
print("The calculation uses the one-loop epsilon expansion around the upper critical dimension (d_c = 4).")
print("-" * 30)
print(f"The scaling relation is: α ≈ (4 - N) / (N + 8) * (4 - d)")
print("\nSubstituting the values:")
print(f"α ≈ ({4} - {N}) / ({N} + {8}) * ({d_c} - {d})")
print(f"α ≈ {numerator} / {denominator} * {epsilon}")
print(f"α ≈ {alpha}")
<<<0.3333333333333333>>>