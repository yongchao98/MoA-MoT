# The task is to find the number of involutions in the group PSU(4, 997).
# An involution is an element g such that g^2=1 and g is not the identity.
# The number of involutions is the sum of the sizes of the conjugacy classes of involutions.
# For PSU(4, q) where q is congruent to 1 mod 4 (997 = 4*249 + 1),
# there are three such classes.
# Summing their sizes leads to the formula:
# N = q^4 * (q+1)^2 * (q^2 - q + 1) * (q^2 + 2*q - 1)

# Set the value for q
q = 997

# Calculate each term in the formula
term1_val = q**4
term2_val = (q + 1)**2
term3_val = q**2 - q + 1
term4_val = q**2 + 2 * q - 1

# Calculate the final result
# The result is an extremely large number, so we will represent it as a product of its terms.
result = term1_val * term2_val * term3_val * term4_val

# Print the equation with the substituted values
print(f"The number of involutions in PSU(4, {q}) is given by the formula:")
print("N = q^4 * (q+1)^2 * (q^2 - q + 1) * (q^2 + 2*q - 1)")
print("\nFor q = 997, the terms are:")
print(f"q^4 = {q}^4 = {term1_val}")
print(f"(q+1)^2 = {q+1}^2 = {term2_val}")
print(f"(q^2 - q + 1) = {q**2 - q + 1} = {term3_val}")
print(f"(q^2 + 2*q - 1) = {q**2 + 2*q - 1} = {term4_val}")
print("\nThe final equation is:")
print(f"N = {term1_val} * {term2_val} * {term3_val} * {term4_val}")
print("\nThe final result is:")
print(result)
