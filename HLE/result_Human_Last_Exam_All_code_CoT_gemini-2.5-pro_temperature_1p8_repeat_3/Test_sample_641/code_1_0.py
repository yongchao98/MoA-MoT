import math

# The prime for the finite field
q = 997

# Calculate the terms of the formula for the number of involutions.
# The formula is N = q^4 * (q^2 - q + 1) * (q^2 + 1) + 1
q_pow_2 = q**2
q_pow_4 = q**4
term1 = q_pow_2 - q + 1
term2 = q_pow_2 + 1

# Calculate the total number of involutions
num_involutions = q_pow_4 * term1 * term2 + 1

# Print the formula and the calculated values
print(f"The number of involutions in PSU(4, q) for q=997 is given by the formula:")
print(f"Number = q^4 * (q^2 - q + 1) * (q^2 + 1) + 1")
print(f"Substituting q = {q}:")
print(f"Number = {q}^4 * ({q}^2 - {q} + 1) * ({q}^2 + 1) + 1")
print(f"Number = {q_pow_4} * ({q_pow_2} - {q} + 1) * ({q_pow_2} + 1) + 1")
print(f"Number = {q_pow_4} * {term1} * {term2} + 1")
print(f"Number = {q_pow_4 * term1 * term2} + 1")
print(f"Final Result = {num_involutions}")
