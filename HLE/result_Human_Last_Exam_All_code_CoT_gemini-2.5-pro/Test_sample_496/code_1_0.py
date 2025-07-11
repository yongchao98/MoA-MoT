# This script calculates the total rank of the SO(4)-equivariant cohomology
# ring of the complement of a 3-dimensional invariant submanifold X, for degrees up to 100.
# The Poincare series for this ring, A, is given by P_A(t) = t^11 / (1-t^4)^2.
# The rank of A in degree k is the coefficient of t^k in the series expansion.
# P_A(t) = sum_{n=0 to inf} (n+1) * t^(4n+11).
# We need to find the total rank for degrees k <= 100.
# This requires summing the ranks for all non-negative integers n
# such that 4n + 11 <= 100.

# Solve for n:
# 4n <= 100 - 11
# 4n <= 89
# n <= 22.25
# So, n ranges from 0 to 22.

total_rank = 0
limit_n = 22

# We calculate the sum S = (0+1) + (1+1) + ... + (22+1)
# which is S = 1 + 2 + ... + 23

# Build a list of the ranks to be summed.
ranks_to_sum = []
for n in range(limit_n + 1):
    rank = n + 1
    ranks_to_sum.append(rank)
    total_rank += rank

# The problem asks to show the final equation.
# The equation for the total rank is the sum of the individual ranks.
print("The total rank is the sum of (n+1) for n from 0 to 22.")
sum_expression = " + ".join(map(str, ranks_to_sum))
print(f"Total Rank = {sum_expression}")

# Print the final computed value.
print(f"The value of the sum is {total_rank}.")
