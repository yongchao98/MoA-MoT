import fractions

# Plan: Calculate the specific heat exponent alpha using the epsilon expansion.

# 1. Define the parameters for the calculation.
# d: spatial dimension
# d_c: upper critical dimension for the n-vector model universality class
# n: number of components of the order parameter (n=1 for the Ising model)
d = 3
d_c = 4
n = 1

# 2. Calculate epsilon (ϵ), the expansion parameter.
epsilon = d_c - d

# 3. Define the numerator and denominator for the alpha formula to the first order in epsilon.
# Formula: α = (4 - n) / (2 * (n + 8)) * ϵ
numerator = 4 - n
denominator = 2 * (n + 8)

# 4. Calculate the value of alpha.
alpha_val = (numerator / denominator) * epsilon

# 5. Use the fractions module to get the exact fractional representation.
fractional_alpha = fractions.Fraction(numerator, denominator) * epsilon

# 6. Print the results, showing each number in the final equation.
print("Calculating the specific heat exponent α for d=3 using the epsilon expansion (to first order).")
print(f"Assuming the Ising model universality class (n={n}).")
print("-" * 20)
print(f"Spatial dimension, d = {d}")
print(f"Upper critical dimension, d_c = {d_c}")
print(f"Epsilon, ε = d_c - d = {d_c} - {d} = {epsilon}")
print(f"Number of order parameter components, n = {n}")
print("-" * 20)
print("The first-order epsilon expansion formula for α is: (4 - n) / (2 * (n + 8)) * ε")
print("Substituting the values:")
print(f"α = ({4} - {n}) / ({2} * ({n} + {8})) * {epsilon}")
print(f"α = {numerator} / {denominator} * {epsilon}")
print(f"α = {fractional_alpha}")
print(f"As a decimal, α ≈ {alpha_val:.4f}")
