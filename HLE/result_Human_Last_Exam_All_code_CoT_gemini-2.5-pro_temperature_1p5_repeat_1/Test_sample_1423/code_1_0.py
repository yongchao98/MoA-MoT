# The number of distinct digits allowed is at most 5.
# To find the maximum possible length of the integer N, we should use the
# maximum number of distinct digits, which is k=5.
k = 5

# As derived in the explanation, the maximum length of a string that satisfies the
# given condition over an alphabet of size k is 2**k - 1.
max_len = 2**k - 1

# Print the equation and the final result.
# The problem asks to output each number in the final equation.
print("The maximum number of digits is found by the equation:")
print(2, "**", k, "-", 1, "=", max_len)
