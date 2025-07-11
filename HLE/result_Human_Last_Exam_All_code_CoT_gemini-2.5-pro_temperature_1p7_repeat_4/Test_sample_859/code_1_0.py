import math

# The problem is about finding the minimal number of edges to add to a graph G'
# to make it 2-edge-connected. The number of edges required depends on 'd'.
# Based on analyzing the graph properties and constraints, the minimal
# number of edges can be determined by a formula involving d.
#
# Let L be the number of "leaf" components in G's block-cut tree structure.
# The number of edges to add is ceil(L/2).
# A detailed combinatorial analysis suggests that the maximum possible number of leaves (worst-case scenario) is L = d + 2.
#
# Therefore, the number of edges required is:
# ceil((d + 2) / 2)
#
# Since d is an even number, let d = 2k for some integer k >= 1.
# The expression becomes:
# ceil((2k + 2) / 2) = ceil(k + 1) = k + 1
#
# Substituting back k = d / 2, the number of edges is:
# (d / 2) + 1
#
# The final code will output this formula. As per the instructions, we must output
# the numbers that appear in the equation. These numbers are 2 and 1.

numerator_of_d = 1 # Not explicitly shown but d is d/1
denominator_of_d = 2
addend = 1

# The problem states d is an even number.
# The minimal number of edges can be expressed as a formula in terms of d.
# The formula is (d/2) + 1.
# We are asked to output each number in the final equation.

print("The minimal number of new edges is given by the equation: d / {} + {}".format(denominator_of_d, addend))

# As a verification for d=2 (the smallest possible even d):
# Edges = 2/2 + 1 = 2.
# For d=4:
# Edges = 4/2 + 1 = 3.
# The result is an integer since d is even.