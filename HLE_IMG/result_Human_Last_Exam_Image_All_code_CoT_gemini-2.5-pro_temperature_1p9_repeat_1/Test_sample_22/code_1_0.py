# This program identifies the indices of knots that are isotopic to the unknot
# from the provided image.

# Based on knot theory, a knot is isotopic to the unknot if it can be
# deformed into a simple circle without cutting or passing it through itself.

# Analysis of each knot:
# K1: 5_1 knot (cinquefoil), not unknot.
# K2: A complex diagram of the unknot. Isotopic to the unknot.
# K3: 3_1 knot (trefoil), not unknot.
# K4: 4_1 knot (figure-eight), not unknot.
# K5: Another complex diagram of the unknot. Isotopic to the unknot.
# K6: 6_1 knot (stevedore), not unknot.

# The indices of the knots isotopic to the unknot are 2 and 5.
unknot_indices = [2, 5]

# Print the list of indices.
print(unknot_indices)