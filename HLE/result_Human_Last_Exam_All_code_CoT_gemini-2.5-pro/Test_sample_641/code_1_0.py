# The user wants to find the number of involutions in PSU(4, 997).
# This corresponds to the projective special unitary group PSU(n, q)
# with n=4 and q=997.

# The formula for the number of involutions in PSU(4, q^2) is:
# I = q^4 * (q^2 + 1) * (q^2 - q + 1)

# Set the value of q
q = 997

# Calculate each term in the formula
term1 = q**4
term2 = q**2 + 1
term3 = q**2 - q + 1

# Calculate the final result by multiplying the terms
# Python's integers handle arbitrary size, so overflow is not an issue.
result = term1 * term2 * term3

# Print the final result in a descriptive equation format as requested.
print("The number of involutions in PSU(4, 997) is calculated using the formula:")
print("Number = q^4 * (q^2 + 1) * (q^2 - q + 1)")
print("\nFor q = 997, the values of the terms are:")
print(f"q^4 = {term1}")
print(f"q^2 + 1 = {term2}")
print(f"q^2 - q + 1 = {term3}")
print("\nThe final equation is:")
print(f"Number = {term1} * {term2} * {term3}")
print(f"\nTotal number of involutions: {result}")
