# The closed-form for the sum was derived using the method of generating functions
# combined with symmetry properties of the summation.
# The final expression obtained for the sum S_n is:
# S_n = (1/8) * 4^n * (n+1)^2 * (63*n^3 + 119*n^2 + 54*n + 8)

# The following code prints this formula and its constituent numerical parts.

print("The closed form for the sum is:")
formula = "S_n = (1/8) * 4**n * (n+1)**2 * (63*n**3 + 119*n**2 + 54*n + 8)"
print(formula)

print("\nAs requested, here are the individual numbers in the final equation:")

fraction_numerator = 1
fraction_denominator = 8
exponential_base = 4
n_plus_1_term_power = 2
poly_coeffs = [63, 119, 54, 8]

print(f"Numerator of the leading fraction: {fraction_numerator}")
print(f"Denominator of the leading fraction: {fraction_denominator}")
print(f"Base of the exponential term: {exponential_base}")
print(f"Exponent of the (n+1) term: {n_plus_1_term_power}")
print(f"Coefficients of the cubic polynomial (from n^3 to the constant term): {poly_coeffs[0]}, {poly_coeffs[1]}, {poly_coeffs[2]}, {poly_coeffs[3]}")