# My thought process to arrive at the values:
# 1. The initial object A is the scale id: Z -> Z. The group G_A is Z.
# 2. The terminal object B in the category of all non-decreasing homomorphisms Z->G is the map to the trivial group, B: Z -> {0}. The group G_B is {0}. Although the problem specifies scales must be "nontrivial", for B to exist, we must assume it is the terminal object in this broader category.
# 3. The given scale S is the inclusion of the integers into the hyperreals, S: Z -> *R. The group G_S is *R.

# Calculation for |S/A|:
# The quotient S/A is formed by the groups G_S / G_A, which is *R / Z.
# The cardinality of this quotient is the cardinality of *R, since Z is countable.
# The cardinality of the hyperreals *R (constructed via an ultrapower over N) is the same as the cardinality of the continuum R.
# |R| = 2^aleph_0 = beth_1.
# So, |S/A| = Beth_1.
card_S_div_A = "Beth_1"

# Calculation for |B/S|:
# The quotient B/S is formed by G_B / Im(h), where h is the map from G_S -> G_B.
# This corresponds to {0} / Im(h: *R -> {0}).
# The only homomorphism from *R to {0} is the zero map, whose image is {0}.
# So the quotient is {0}/{0}, which is the trivial group.
# The cardinality of the trivial group is 1.
card_B_div_S = 1

# Calculation for H_1(B/A, Q):
# The space B/A is formed by G_B / Im(h), where h is the map from G_A -> G_B.
# This corresponds to {0} / Im(h: Z -> {0}).
# The map h is the zero map, and its image is {0}.
# The space B/A is {0}/{0}, which is a single point space.
# The first homology group of a point, H_1({*}, Q), is the trivial group {0}.
# The cardinality of the trivial group is 1.
card_H1_B_div_A = 1

# The problem asks for the answer in the format of space-separated values.
print(f"{card_S_div_A} {card_B_div_S} {card_H1_B_div_A}")
