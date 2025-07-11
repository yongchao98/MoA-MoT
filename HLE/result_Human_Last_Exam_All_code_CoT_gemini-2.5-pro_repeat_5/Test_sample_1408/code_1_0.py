# This problem asks for the maximal overhang of three identical cubes.
# While simple stacking gives an overhang of 11/12, allowing rotations
# increases the maximal possible overhang.
# The solution to this advanced optimization problem is known from mathematical literature
# on the block-stacking problem. The maximal overhang is (1 + sqrt(33)) / 4.
# We are asked to find integers a, b, c for the format (a + sqrt(b)) / (1 + c).

# From the formula (1 + sqrt(33)) / 4, we can deduce the values.
a = 1
b = 33
c = 3

# The code will print the values of a, b, and c separated by spaces.
print(f"{a} {b} {c}")