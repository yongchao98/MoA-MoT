import math

# This script solves the Theatre Square problem and calculates the answer for the specific inputs from Question 4.

# Input values from Question 4
n = 4 * 10**9
m = 4 * 10**9
a = 1

# To find the number of flagstones needed for one dimension, we calculate the ceiling of the division.
# Using integer arithmetic, ceiling(x/y) can be calculated as (x + y - 1) // y.
# This method avoids using if-statements or floating-point math.

# Number of flagstones along the n-side
na = (n + a - 1) // a

# Number of flagstones along the m-side
ma = (m + a - 1) // a

# Total number of flagstones is the product
total_flagstones = na * ma

# As requested, here are the numbers in the final equation:
print(f"To solve for n={int(n)}, m={int(m)}, a={int(a)}:")
print(f"Flagstones for n-side (na): {na}")
print(f"Flagstones for m-side (ma): {ma}")
print(f"Final Equation: {na} * {ma} = {total_flagstones}")