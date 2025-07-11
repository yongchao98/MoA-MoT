import fractions

# The problem is to calculate the sum of 1/2^(i+j) over a set S.
# Our derivation shows that the set S consists of three disjoint sets of pairs (i, j):
# 1. j = i, for all positive integers i.
# 2. j = 2i, for all positive integers i.
# 3. i = 2j, for all positive integers j.

# We calculate the sum for each case. Each sum is an infinite geometric series.
# The sum of a geometric series sum_{n=1 to inf} ar^(n-1) is a/(1-r).
# For sum_{n=1 to inf} r^n, the first term is r, so the sum is r/(1-r).

# Case 1: j = i
# The sum is sum_{i=1 to inf} 1/2^(i+i) = sum_{i=1 to inf} (1/4)^i.
# The common ratio is r = 1/4.
r1 = fractions.Fraction(1, 4)
sum1 = r1 / (1 - r1)

# Case 2: j = 2i
# The sum is sum_{i=1 to inf} 1/2^(i+2i) = sum_{i=1 to inf} (1/8)^i.
# The common ratio is r = 1/8.
r2 = fractions.Fraction(1, 8)
sum2 = r2 / (1 - r2)

# Case 3: i = 2j
# The sum is sum_{j=1 to inf} 1/2^(2j+j) = sum_{j=1 to inf} (1/8)^j.
# This is identical to Case 2.
sum3 = sum2

# The total sum is the sum of the results from the three cases.
total_sum = sum1 + sum2 + sum3

# As requested, here is the final equation with each number printed.
print(f"{sum1} + {sum2} + {sum3} = {total_sum}")
