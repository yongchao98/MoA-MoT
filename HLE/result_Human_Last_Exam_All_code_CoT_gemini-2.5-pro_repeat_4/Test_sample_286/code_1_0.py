import math

# The problem reduces to finding the smallest integer M satisfying the inequality:
# M >= 50000 * log10(2)

# Define the coefficient in the inequality
coefficient = 50000

# Calculate log10(2)
log10_2 = math.log10(2)

# Calculate the lower bound for M
m_lower_bound = coefficient * log10_2

# Since M must be an integer, we take the ceiling of this value.
M = math.ceil(m_lower_bound)

# Output the components of the final calculation and the result
print("The problem simplifies to the inequality: M >= N * log10(2)")
print(f"Here, the effective number of steps N is {coefficient}.")
print(f"The value of log10(2) is approximately: {log10_2}")
print(f"Thus, M must be greater than or equal to {coefficient} * {log10_2} = {m_lower_bound}")
print(f"The smallest positive integer M satisfying this is the ceiling of {m_lower_bound}, which is {M}.")
