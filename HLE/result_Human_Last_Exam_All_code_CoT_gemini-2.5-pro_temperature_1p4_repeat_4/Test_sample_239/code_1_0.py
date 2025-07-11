import math

#
# Computer:
# This script solves the Theatre Square paving problem.
#
# Problem: Find the minimum number of square flagstones of size a*a
# needed to pave a rectangular area of size n*m.
#

# Input values from the problem description's final question.
n = 4 * 10**9
m = 4 * 10**9
a = 1

# Plan:
# To cover a dimension of length 'L' with tiles of side 's', we need
# ceiling(L/s) tiles. The ceiling function rounds a number up to the
# nearest integer. For example, ceiling(10/3) = ceiling(3.33) = 4.
#
# In integer arithmetic, the ceiling of a division (e.g., n / a) can be
# calculated efficiently without using floating-point math by the formula:
# (n + a - 1) // a
#
# We apply this formula to both the n and m dimensions.

# Calculate the number of flagstones needed along the n-dimension.
flagstones_for_n = (n + a - 1) // a

# Calculate the number of flagstones needed along the m-dimension.
flagstones_for_m = (m + a - 1) // a

# The total number of flagstones is the product of the amounts for each dimension.
total_flagstones = flagstones_for_n * flagstones_for_m

# Print the final equation as requested.
print(f"{flagstones_for_n} * {flagstones_for_m} = {total_flagstones}")