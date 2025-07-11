# The exponent of N in the threshold for the sum's magnitude
level_set_exponent_numerator = 3
level_set_exponent_denominator = 8

# The power 'q' used in the L^q norm and Chebyshev's inequality.
# We choose q=4 based on a known harmonic analysis theorem.
q = 4

# According to the theorem, the L^4 norm of the maximal function is bounded
# by N^epsilon, so its contribution to the exponent is 0 in the limit.
integral_bound_exponent = 0

# Calculate alpha using the formula: alpha = integral_exponent - level_set_exponent * q
alpha = integral_bound_exponent - (level_set_exponent_numerator / level_set_exponent_denominator) * q

# Print out the equation
print(f"The equation for alpha is:")
print(f"alpha = {integral_bound_exponent} - ({level_set_exponent_numerator}/{level_set_exponent_denominator}) * {q}")
print(f"alpha = {alpha}")
