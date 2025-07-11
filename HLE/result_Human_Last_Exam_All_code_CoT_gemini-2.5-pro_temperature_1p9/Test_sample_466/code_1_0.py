# Set the parameters for the finite field and vector space dimension
q = 11
n = 3

# Calculate the components of the formula for the order of GL(n, q)
# |GL(n,q)| = (q^n - 1)(q^n - q)...(q^n - q^(n-1))
term1 = q**n - 1
term2 = q**n - q
term3 = q**n - q**2

# Calculate the final result
result = term1 * term2 * term3

# Print the final equation with all numbers, as requested.
# The number of adjunctions is equal to the number of invertible 3x3 matrices over F_11,
# which is the order of the group GL(3, 11).
print(f"The number of adjunctions is given by the equation:")
print(f"({11**3} - 1) * ({11**3} - 11) * ({11**3} - {11**2}) = {term1} * {term2} * {term3} = {result}")
