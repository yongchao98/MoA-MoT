import sys
# This is a reasoning problem based on knot theory.
# The analysis is performed by visually inspecting the knot diagrams
# and identifying which ones can be simplified to the unknot using Reidemeister moves.

# K_1 is the 6_1 knot (Stevedore knot), which is not the unknot.
# K_2 can be untangled using Reidemeister moves (one type I, two type II). It is the unknot.
# K_3 is the 3_1 knot (trefoil knot), which is not the unknot.
# K_4 is the 4_1 knot (figure-eight knot), which is not the unknot.
# K_5 can be untangled using Reidemeister moves (one type I, two type II). It is the unknot.
# K_6 is the 5_1 knot (cinquefoil knot), which is not the unknot.

# The indices of the knots isotopic to the unknot are 2 and 5.
unknot_indices = [2, 5]

# The problem asks to return a list of indices.
# We will print the list to the standard output.
print(unknot_indices)
