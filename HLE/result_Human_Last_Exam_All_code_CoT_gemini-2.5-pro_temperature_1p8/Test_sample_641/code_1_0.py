# The value of q for PSU(4, q)
q = 997

# Calculate the terms of the formula
q_pow_4 = q**4
q_sq_plus_1 = q**2 + 1
q_sq_minus_q_plus_1 = q**2 - q + 1

# Calculate the numerator of the formula, which is the number of
# non-central involutions in SU(4, q)
numerator = q_pow_4 * q_sq_plus_1 * q_sq_minus_q_plus_1

# The number of involutions in PSU(4, q) is the numerator divided by 2
num_involutions = numerator // 2

# Print the final equation with the computed result.
print(f"Number of involutions = ({q}^4 * ({q}^2 + 1) * ({q}^2 - {q} + 1)) / 2")
print(f"= ({q_pow_4} * {q_sq_plus_1} * {q_sq_minus_q_plus_1}) / 2")
print(f"= {num_involutions}")