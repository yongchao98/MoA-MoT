# The dimension of the ambient space for the curve is d.
# The curve is in R^3, so d=3.
d = 3

# The formula for the sharp l^2 decoupling exponent (the critical exponent p_c)
# for a non-degenerate curve in R^d is:
# p_c = d * (d + 1) / (d - 1)
# This result is a consequence of the l^2 decoupling theorem by Bourgain, Demeter, and Guth.

# We calculate the numerator and denominator of the formula.
numerator = d * (d + 1)
denominator = d - 1

# Calculate the final value for the exponent.
p_c = numerator / denominator

# Print the final equation with all the numbers, as requested.
print(f"The sharp l^2 decoupling exponent p_c is calculated as:")
print(f"p_c = {d} * ({d} + 1) / ({d} - 1)")
print(f"p_c = {numerator} / {denominator}")
print(f"p_c = {p_c}")
