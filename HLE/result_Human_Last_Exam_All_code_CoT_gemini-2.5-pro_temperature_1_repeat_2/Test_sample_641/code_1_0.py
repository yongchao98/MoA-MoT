# Plan:
# 1. Define the value of q.
# 2. Calculate the components of the formula: q^4, (q^2+1), and (q^2-q+1).
# 3. Multiply these components together.
# 4. Divide by 2 to get the final count.
# 5. Print the result in a formatted equation.

# Set the value for q
q = 997

# Calculate the terms of the formula
# Note: Python's integers handle arbitrarily large numbers, so overflow is not an issue.
term1 = q**4
term2 = q**2 + 1
term3 = q**2 - q + 1

# Calculate the number of non-central involutions in SU(4, q)
numerator = term1 * term2 * term3

# The number of involutions in PSU(4, q) is half of this value.
# The result is guaranteed to be an integer because q is odd, making q^2+1 even.
num_involutions = numerator // 2

# Print the final equation with all the numbers
print(f"({q}**4 * ({q}**2 + 1) * ({q}**2 - {q} + 1)) / 2 =")
print(f"({term1} * {term2} * {term3}) / 2 =")
print(f"{numerator} / 2 =")
print(num_involutions)
