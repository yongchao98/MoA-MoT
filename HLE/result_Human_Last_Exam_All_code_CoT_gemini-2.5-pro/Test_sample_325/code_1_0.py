# The dimension of the ambient space for the curve is d.
d = 3

# The formula for the sharp l^2 decoupling exponent p for a non-degenerate
# curve in R^d is: p = 2 * (d + 1) / (d - 1).

# Numerator of the formula
numerator = 2 * (d + 1)

# Denominator of the formula
denominator = d - 1

# Calculate the final exponent p
p = numerator / denominator

# Print the step-by-step calculation of the exponent
print("The formula for the sharp l^2 decoupling exponent p is:")
print("p = 2 * (d + 1) / (d - 1)")
print("")
print(f"The given curve is in R^3, so the dimension d = {d}.")
print("Substituting d = 3 into the formula:")
print(f"p = 2 * ({d} + 1) / ({d} - 1)")
print(f"p = 2 * {d + 1} / {d - 1}")
print(f"p = {numerator} / {denominator}")
print(f"p = {p}")