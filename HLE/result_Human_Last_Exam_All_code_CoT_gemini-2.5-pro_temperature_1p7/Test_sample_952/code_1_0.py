# The size of the square table.
N = 2024

# The problem is to find the largest k such that any arrangement of k diamonds is not "stuck".
# This is equivalent to k_min - 1, where k_min is the minimum number of diamonds
# needed to form a "stuck" configuration.

# For an N x N grid where N is even, the minimum size of a stuck configuration
# is achieved by placing diamonds on every other row and column.
# The number of such diamonds is (N/2) * (N/2).
half_N = N // 2
k_min = half_N * half_N

# The largest number of diamonds for which every arrangement is guaranteed to be "mobile" is k_min - 1.
result = k_min - 1

# Print the calculation and the final answer.
print(f"The calculation is based on the size of the minimal 'stuck' configuration, which is ({N}/2)^2.")
print(f"So, the largest value k is ({N}/2)^2 - 1.")
print(f"{half_N} * {half_N} - 1 = {result}")
