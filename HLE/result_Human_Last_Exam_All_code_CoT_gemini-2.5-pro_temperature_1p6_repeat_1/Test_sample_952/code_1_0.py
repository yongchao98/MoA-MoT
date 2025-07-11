# The size of the square table.
N = 2024

# The problem is to find the largest k such that any arrangement of k diamonds
# has a movable piece. This is equivalent to M-1, where M is the size of the
# smallest "frozen" arrangement (an arrangement with no movable pieces).

# A minimal frozen arrangement can be constructed by placing diamonds
# on cells (i, j) where both i and j are even.
# For an N x N grid (where N is even), the number of even indices (0, 2, ..., N-2)
# is N / 2.
side_count = N // 2

# The total number of diamonds in this minimal frozen arrangement (M) is:
M = side_count * side_count

# The largest k for which a move is always possible is M - 1.
k = M - 1

# The problem asks to output each number in the final equation.
print(f"The size of the minimal frozen arrangement is M = ({N} / 2)^2 = {side_count}^2 = {M}.")
print(f"The largest value k is M - 1.")
print(f"k = {M} - 1 = {k}")
