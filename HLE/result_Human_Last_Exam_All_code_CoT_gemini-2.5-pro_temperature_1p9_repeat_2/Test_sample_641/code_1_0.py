# Set the value for q from PSU(4, q)
q = 997

# The formula for the number of involutions in PSU(4, q) for q = 1 (mod 4) is:
# 1 + N2 / 2, where N2 is the number of non-central involutions of a specific type in SU(4, q).
# N2 = q^4 * (q^2 + 1) * (q^2 - q + 1)

# Perform the calculation using Python's support for large integers.
q_squared = q**2
q_fourth = q**4

# These are the terms in the formula for N2
term1 = q_squared + 1
term2 = q_squared - q + 1

# Calculate the number of non-central involutions in SU(4, q)
N2 = q_fourth * term1 * term2

# The number of involutions in PSU(4, q) is 1 (the central one) + N2 / 2
num_involutions = 1 + N2 // 2

# Print the final equation and the result as requested.
print("The number of involutions in PSU(4, 997) is calculated by the formula:")
print("Result = 1 + (q^4 * (q^2 + 1) * (q^2 - q + 1)) / 2")
print("\nSubstituting q = 997:")
print(f"Result = 1 + ({q}**4 * ({q}**2 + 1) * ({q}**2 - {q} + 1)) // 2")
print("\nCalculating the values:")
print(f"Result = 1 + ({q_fourth} * {term1} * {term2}) // 2")
print("\nFinal Answer:")
print(f"Number of involutions = {num_involutions}")
