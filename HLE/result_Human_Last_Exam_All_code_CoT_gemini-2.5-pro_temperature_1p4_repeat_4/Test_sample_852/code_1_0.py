# Based on the mathematical theorem by Kedlaya, the smallest size of a finite
# Abelian group G containing a maximal sum-free set S with |k(S)| > 2|S|
# is determined by specific properties of G.
# The theorem states this is possible if and only if the subgroup K of elements
# of order at most 2 has size |K| >= 8, and the quotient group G/K has
# exponent 3.
# This implies |G| must be a multiple of |K| * |G/K|, which in turn must be
# a multiple of 8 * 3 = 24.
# The smallest such group is Z₂ x Z₂ x Z₂ x Z₃, with order 24.

smallest_size = 24
print(smallest_size)