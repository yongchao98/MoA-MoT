# The number of distinct digits allowed.
k = 5

# The maximum length of a valid string with k distinct digits
# follows the formula f(k) = 2^k - 1.
# We calculate this for k=5.

# Calculate 2 to the power of k
power_of_2 = 2**k

# The result is one less than that
result = power_of_2 - 1

# Print the final equation and the result, as requested.
print(f"The maximum possible number of digits is given by the formula 2^k - 1.")
print(f"For k = {k}:")
print(f"2^{k} - 1 = {power_of_2} - 1 = {result}")
