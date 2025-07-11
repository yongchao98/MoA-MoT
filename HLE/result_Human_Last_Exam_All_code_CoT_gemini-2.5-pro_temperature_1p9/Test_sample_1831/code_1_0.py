# A program to calculate the number of equivalence classes of endomorphisms on a 4-element set.

# c(n) is the number of non-isomorphic connected functional graphs on n vertices.
# Based on enumeration of structures (cycle with rooted trees attached).
c1 = 1  # A single fixed point.
c2 = 2  # A 2-cycle, or a path of length 2 into a fixed point.
c3 = 4  # A 3-cycle, a 2-cycle with one leaf, or two types of 3-vertex rooted trees.
c4 = 9  # A 4-cycle, a 3-cycle with one leaf, 3 types on a 2-cycle, or 4 types of 4-vertex rooted trees.

# We sum the number of possibilities for each partition of the integer 4.

# Partition [4]: one component of size 4.
# The number of structures is c(4).
num_part_4 = c4
print(f"Number of structures for partition [4]: {num_part_4}")

# Partition [3, 1]: one component of size 3 and one of size 1.
# The number of structures is c(3) * c(1).
num_part_3_1 = c3 * c1
print(f"Number of structures for partition [3, 1]: {num_part_3_1}")

# Partition [2, 2]: two components of size 2.
# We are choosing 2 items from c(2) types, with replacement.
# This is given by the multiset coefficient formula binom(n+k-1, k), with n=c(2) and k=2.
# The types can be (type1, type1), (type2, type2), (type1, type2).
num_part_2_2 = 3
print(f"Number of structures for partition [2, 2]: {num_part_2_2}")

# Partition [2, 1, 1]: one component of size 2, two of size 1.
# We choose one structure of size 2 (c(2) ways) and the rest is fixed (one way for two size-1 components).
num_part_2_1_1 = c2
print(f"Number of structures for partition [2, 1, 1]: {num_part_2_1_1}")

# Partition [1, 1, 1, 1]: four components of size 1.
# All components are fixed points. Only one way.
num_part_1_1_1_1 = 1
print(f"Number of structures for partition [1, 1, 1, 1]: {num_part_1_1_1_1}")

# The total number of equivalence classes is the sum of these possibilities.
total = num_part_4 + num_part_3_1 + num_part_2_2 + num_part_2_1_1 + num_part_1_1_1_1
print(f"\nTotal number of equivalence classes is the sum:")
print(f"{num_part_4} + {num_part_3_1} + {num_part_2_2} + {num_part_2_1_1} + {num_part_1_1_1_1} = {total}")
