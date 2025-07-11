# This script calculates the total rank of the SO(4)-equivariant cohomology
# ring of SO(4)\X up to degree 100, where X is a 3-dimensional invariant
# submanifold identified as SO(3).

# The Poincaré series for the cohomology ring A is P_A(t) = (t^3 + t^6) / (1-t^4)^2.
# We expand this as a power series:
# P_A(t) = (t^3 + t^6) * sum_{j=0 to inf} (j+1)*t^(4j)
#        = sum_{j=0 to inf} (j+1)*t^(4j+3) + sum_{j=0 to inf} (j+1)*t^(4j+6)
#
# We need to sum the coefficients for degrees up to 100.

# For the first sum, we need degree 4j+3 <= 100, which means 4j <= 97, so j <= 24.
# The coefficients are (j+1) for j from 0 to 24.
# This is the sum of integers from 1 to 25.
max_j1 = (100 - 3) // 4
sum1 = sum(j + 1 for j in range(max_j1 + 1))

# For the second sum, we need degree 4j+6 <= 100, which means 4j <= 94, so j <= 23.
# The coefficients are (j+1) for j from 0 to 23.
# This is the sum of integers from 1 to 24.
max_j2 = (100 - 6) // 4
sum2 = sum(j + 1 for j in range(max_j2 + 1))

# The total rank is the sum of these two parts.
total_rank = sum1 + sum2

print(f"The calculation is based on the sum of two series.")
print(f"Sum of ranks from the first series (t^3 term): {sum1}")
print(f"Sum of ranks from the second series (t^6 term): {sum2}")
print(f"The total rank is the sum of these two values.")
print(f"Total Rank = {sum1} + {sum2} = {total_rank}")