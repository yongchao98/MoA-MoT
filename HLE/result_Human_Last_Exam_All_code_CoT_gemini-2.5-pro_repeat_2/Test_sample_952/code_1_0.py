# The side length of the square table.
N = 2024

# The problem is to find the largest k such that any arrangement of k diamonds has a movable diamond.
# This is equivalent to finding M-1, where M is the minimum number of diamonds in an "all-stuck"
# arrangement (an arrangement where no diamond can be moved).
#
# A minimal all-stuck arrangement can be constructed by placing diamonds on cells (i, j)
# where both i and j are even. For an N x N grid where N is even, this means placing diamonds
# on a sub-grid of size (N/2) x (N/2).
#
# The number of diamonds in this minimal all-stuck arrangement is M = (N/2)^2.
# This is known to be the minimum size for such a configuration.
# The value k we seek is M - 1.

# Calculation steps:
# 1. Get the side length of the sub-grid.
half_n = N // 2

# 2. Calculate the number of diamonds in the minimal all-stuck configuration (M).
min_stuck_diamonds = half_n ** 2

# 3. The result k is M - 1.
k = min_stuck_diamonds - 1

# Output the reasoning and the final equation with all numbers.
# The explanation is included in the printout for clarity.
print(
    f"The problem is to find the largest k such that any arrangement of k diamonds has at least one movable diamond.\n"
    f"This value is one less than the minimum number of diamonds (M) in an 'all-stuck' configuration.\n"
    f"For a grid of size N = {N}, the minimum size of an 'all-stuck' configuration is M = (N/2)^2.\n"
    f"Calculation of M: ({N} / 2)^2 = {half_n}^2 = {min_stuck_diamonds}.\n"
    f"The largest value k is M - 1.\n"
    f"Final Equation: ({N} / 2)^2 - 1 = {k}"
)
