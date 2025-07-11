import sys
# The user wants a list of indices i for which the knot K_i is the unknot.
# Based on visual inspection and knowledge of knot theory:
# K_1 is the stevedore knot (6_1), which is not the unknot.
# K_2 is a complex representation of the unknot.
# K_3 simplifies to the unknot via one Reidemeister II move.
# K_4 is the figure-eight knot (4_1), which is not the unknot.
# K_5 simplifies to the trefoil knot, so it is not the unknot.
# K_6 is the trefoil knot (3_1), which is not the unknot.
# Therefore, only knots 2 and 3 are the unknot.

unknot_indices = [2, 3]

# Print the final list of indices.
# The list [2, 3] represents that knots K_2 and K_3 are isotopic to the unknot.
print(unknot_indices)