# Parameters for the calculation
n = 1      # Number of spin components for the Ising model
d = 3      # Spatial dimension
d_c = 4    # Upper critical dimension

# Calculate epsilon
epsilon = d_c - d

# The formula for the specific heat exponent alpha to first order in epsilon is:
# alpha = (4 - n) / (2 * (n + 8)) * epsilon

# Perform the calculation
numerator = 4 - n
denominator_factor_1 = 2
denominator_factor_2 = n + 8
alpha = (numerator / (denominator_factor_1 * denominator_factor_2)) * epsilon

print("Calculation of the specific heat exponent α:")
print(f"Using n={n} (Ising model), d={d}, d_c={d_c}, which gives epsilon = {epsilon}")
print("\nFormula: α = (4 - n) * ε / (2 * (n + 8))")
print("Substituting the values:")
print(f"α = (4 - {n}) * {epsilon} / (2 * ({n} + 8))")
print(f"α = ({numerator}) * {epsilon} / (2 * ({denominator_factor_2}))")
print(f"α = {numerator} / {denominator_factor_1 * denominator_factor_2}")
print(f"α ≈ {alpha:.4f}")
