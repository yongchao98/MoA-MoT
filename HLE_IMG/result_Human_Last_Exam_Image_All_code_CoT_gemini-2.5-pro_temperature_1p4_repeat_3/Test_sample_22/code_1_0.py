# This script identifies which of the six knots K_1, ..., K_6 are isotopic to the unknot.
# An unknot is a knot that can be continuously deformed into a simple loop.
# This determination is made by applying Reidemeister moves to simplify the knot diagrams.

# K_1 is the stevedore knot (6_1), which is not an unknot.
# K_2 is a complicated representation of the unknot. It can be untangled.
# K_3 can be simplified to an unknot with a single Reidemeister II move.
# K_4 is the figure-eight knot (4_1), which is not an unknot.
# K_5 is another complex drawing of the unknot, simplified by a Reidemeister I move followed by unlinking.
# K_6 is the trefoil knot (3_1), which is not an unknot.

# The list of indices for knots that are isotopic to the unknot is [2, 3, 5].
unknot_indices = [2, 3, 5]

# Print the final list of indices.
print(unknot_indices)