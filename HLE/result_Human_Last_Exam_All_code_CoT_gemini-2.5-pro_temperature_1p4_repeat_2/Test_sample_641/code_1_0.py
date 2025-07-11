# Set the value for q as specified in the problem
q = 997

# The formula for the number of involutions in PSU(4, q) for q = 1 (mod 4) is:
# N = q^4 * (q+1)^2 * (q^2 - q + 1) * (q^2 + 1)

# Calculate each term of the formula
term1 = q**4
term2 = (q + 1)**2
term3 = q**2 - q + 1
term4 = q**2 + 1

# Calculate the final result by multiplying the terms
result = term1 * term2 * term3 * term4

# Print the calculation steps and the final answer
print(f"The number of involutions in PSU(4, 997) is calculated using the formula:")
print(f"N = q^4 * (q+1)^2 * (q^2 - q + 1) * (q^2 + 1)")
print(f"\nFor q = 997, the terms are:")
print(f"q^4 = {q}^4 = {term1}")
print(f"(q+1)^2 = {q+1}^2 = {term2}")
print(f"q^2 - q + 1 = {q}^2 - {q} + 1 = {term3}")
print(f"q^2 + 1 = {q}^2 + 1 = {term4}")
print(f"\nThe final equation is:")
print(f"{term1} * {term2} * {term3} * {term4} = {result}")