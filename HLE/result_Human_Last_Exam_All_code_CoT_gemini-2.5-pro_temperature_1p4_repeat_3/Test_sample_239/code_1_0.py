import math

# Computer:

# For this demonstration, we will use the specific inputs from Question 4.
# In a general case, these values would be read from user input, for example:
# n, m, a = map(int, input().split())
n = 4000000000
m = 4000000000
a = 1

# To find the number of flagstones needed along one dimension, we need to calculate
# the ceiling of the dimension's length divided by the flagstone's side length.
# In integer arithmetic, ceil(x / y) can be calculated as (x + y - 1) // y.
# This avoids using floating-point math and is more efficient.

# Calculate flagstones needed for the 'n' dimension.
flagstones_n = (n + a - 1) // a

# Calculate flagstones needed for the 'm' dimension.
flagstones_m = (m + a - 1) // a

# The total number of flagstones is the product of the two dimensions.
total_flagstones = flagstones_n * flagstones_m

# As requested, we print the full equation showing each number.
# The 'f-string' formatting is used to embed the variable values directly into the string.
print(f"ceil({n} / {a}) * ceil({m} / {a}) = {flagstones_n} * {flagstones_m} = {total_flagstones}")