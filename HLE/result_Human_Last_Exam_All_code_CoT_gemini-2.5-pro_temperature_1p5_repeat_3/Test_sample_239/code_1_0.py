import math

# --- Inputs from Question 4 ---
n = 4 * 10**9
m = 4 * 10**9
a = 1

# --- Calculation ---
# To find the number of flagstones needed for one dimension, we need to
# calculate the ceiling of the division. For positive integers,
# ceil(x / y) can be calculated using integer division as (x + y - 1) // y.

# Number of flagstones for the n-side
na = (n + a - 1) // a

# Number of flagstones for the m-side
ma = (m + a - 1) // a

# Total number of flagstones
total_flagstones = na * ma

# --- Output ---
# As requested, here is the final equation with each number.
print(f"To solve the Theatre Square problem for n={n}, m={m}, a={a}:")
print(f"Number of stones needed = ({n} + {a} - 1) // {a} * ({m} + {a} - 1) // {m}")
print(f"Equation: {na} * {ma} = {total_flagstones}")