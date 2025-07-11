# The reasoning above establishes a lower bound and an upper bound for the cardinality.

# 1. A theorem by Waraszkiewicz proves that any decomposable continuum
#    must have at least two regular proper subcontinua.
lower_bound = 2

# 2. An example construction (the union of two pseudo-arcs at a single point)
#    provides a decomposable continuum that has exactly two regular
#    proper subcontinua.
achievable_cardinality = 2

# Since the lower bound is achievable, the smallest possible cardinality is 2.
smallest_possible_cardinality = 2

# There is no equation, so we print the final resulting number.
print(smallest_possible_cardinality)