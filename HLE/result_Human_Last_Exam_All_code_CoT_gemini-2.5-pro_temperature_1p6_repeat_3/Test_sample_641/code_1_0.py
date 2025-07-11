# The prime for the finite field in PSU(4, q)
q = 997

# The formula for the number of involutions in PSU(4, q) for q such that gcd(4, q+1)=2 is:
# N = (q**5 * (q**2 - q + 1) * (q**2 + 2*q - 1)) / 2
# We can rewrite this to use integer arithmetic throughout.

print(f"Calculating the number of involutions for q = {q}")
print("The formula is: (q^5 * (q^2 - q + 1) * (q^2 + 2q - 1)) / 2")
print("-" * 30)

# Calculate each term of the formula
term1 = q**5
term2 = q**2 - q + 1
term3_num = q**2 + 2 * q - 1

# Since term3_num is guaranteed to be even, we can perform integer division
term3 = term3_num // 2

print(f"The equation becomes: {term1} * {term2} * ({term3_num} / 2)")
print(f"Which simplifies to: {term1} * {term2} * {term3}")
print("-" * 30)

# Calculate the final result
number_of_involutions = term1 * term2 * term3

print(f"The total number of involutions in PSU(4, {q}) is:")
print(number_of_involutions)