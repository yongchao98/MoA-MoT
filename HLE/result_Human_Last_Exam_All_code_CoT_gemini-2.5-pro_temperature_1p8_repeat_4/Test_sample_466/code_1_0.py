# Calculate the number of internal adjunctions
# This is equal to the order of the general linear group GL(3, F_11)

q = 11
n = 3

# Calculate the terms
term1 = q**n - q**0
term2 = q**n - q**1
term3 = q**n - q**2

# Calculate the total number
result = term1 * term2 * term3

# Print the final equation with all numbers
print(f"Number of adjunctions = ({q}^{n} - {q}^0) * ({q}^{n} - {q}^1) * ({q}^{n} - {q}^2)")
print(f"= ({q**n} - {q**0}) * ({q**n} - {q**1}) * ({q**n} - {q**2})")
print(f"= {term1} * {term2} * {term3}")
print(f"= {result}")