# The problem asks for the minimum total number of edges in a topologically
# nontrivial link with three components on the 3D integer lattice.

# Based on research in knot theory, the minimal such link is a "3-chain,"
# where three knots are linked sequentially (K1 links K2, K2 links K3).
# This configuration can be achieved with three knots of 8 edges each.
# An 8-edge knot can be formed by a 2x2 square on the lattice.

# Length of the first knot (K1)
knot1_length = 8

# Length of the second knot (K2)
knot2_length = 8

# Length of the third knot (K3)
knot3_length = 8

# The total number of edges is the sum of the lengths of the three knots.
total_length = knot1_length + knot2_length + knot3_length

# The final code must output each number in the final equation.
print(f"{knot1_length} + {knot2_length} + {knot3_length} = {total_length}")