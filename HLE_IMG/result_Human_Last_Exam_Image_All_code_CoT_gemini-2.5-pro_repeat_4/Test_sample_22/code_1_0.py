# This script identifies the knots that are isotopic to the unknot from the provided image.
# Based on visual analysis using Reidemeister moves:
# K1 is the Stevedore knot (6_1), not the unknot.
# K2 can be untangled using one Reidemeister II move and one Reidemeister I move. It is the unknot.
# K3 is the square knot, not the unknot.
# K4 is the figure-eight knot (4_1), not the unknot.
# K5 can be untangled using one Reidemeister II move and one Reidemeister I move. It is the unknot.
# K6 is the cinquefoil knot (5_1), not the unknot.

# The list of indices for knots that are isotopic to the unknot.
unknot_indices = [2, 5]

# Print the final list of indices.
print(unknot_indices)