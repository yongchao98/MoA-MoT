import sys

# --- Parameters ---
# The spatial dimension of the system
d = 3
# The upper critical dimension for the O(n) model
dc = 4
# The number of components of the order parameter.
# We assume the Ising model, where n = 1.
# Other common models: XY model (n=2), Heisenberg model (n=3)
n = 1

# --- Calculation ---
# 1. Calculate the expansion parameter, epsilon
epsilon = dc - d

# 2. Calculate the specific heat exponent alpha to the first order in epsilon.
# The general formula for the O(n) model is: alpha = (4 - n) * epsilon / (2 * (n + 8))
# We handle the case where the denominator might be zero, although not for n >= 0.
if (n + 8) == 0:
    print("Error: Division by zero. The number of components 'n' cannot be -8.", file=sys.stderr)
    sys.exit(1)

alpha_numerator = (4 - n) * epsilon
alpha_denominator = 2 * (n + 8)
alpha = alpha_numerator / alpha_denominator

# --- Output ---
print("Calculation of the specific heat scaling exponent α using the epsilon expansion.")
print(f"The calculation is performed for a system with spatial dimension d = {d}.")
print(f"The upper critical dimension for this class of models is dc = {dc}.")
print(f"The expansion parameter is therefore ε = dc - d = {dc} - {d} = {epsilon}.")
print(f"We assume the Ising model, so the number of components of the order parameter is n = {n}.")
print("\nThe general first-order formula for α is: α = (4 - n) * ε / (2 * (n + 8))")
print("\nSubstituting the values into the formula:")
print(f"α = (4 - {n}) * {epsilon} / (2 * ({n} + 8))")
print(f"α = {alpha_numerator} / {alpha_denominator}")
print(f"α ≈ {alpha}")