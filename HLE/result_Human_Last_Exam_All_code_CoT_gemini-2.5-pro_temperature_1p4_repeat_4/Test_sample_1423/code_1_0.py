# The number of distinct digits allowed.
k = 5

# The maximum possible number of digits in N is given by the formula 2^k - 1.
# This formula is derived from the properties of strings where every substring
# must have a character that appears exactly once.
# Let's calculate this value for k=5.
max_digits = 2**k - 1

# Print the final equation and the result.
print(f"2**{k} - 1 = {max_digits}")