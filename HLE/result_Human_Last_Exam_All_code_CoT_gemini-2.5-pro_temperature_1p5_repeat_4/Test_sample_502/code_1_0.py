# The order of the group G.
group_order = 10000

# The dimension of the vector space V the group acts on.
# The polynomial ring is R = C[x_1, ..., x_10].
space_dimension = 10

# According to a theorem from invariant theory, the dimension of the coinvariant
# algebra R/I is equal to the order of the group effectively acting on the space.
# If the action is defined by a representation rho: G -> GL(V), the dimension is |rho(G)|.
# The relationship is dim(R/I) = |G| / |ker(rho)|.

# To maximize this dimension, we must minimize the order of the kernel of the representation.
# The minimum possible order for a kernel is 1, which corresponds to a faithful representation.
min_kernel_order = 1

# For this maximum to be achievable, there must exist a group of order 10000
# that has a faithful representation of dimension 10.
# The abelian group G = (Z/10Z)^4 is such a group.
# It has order 10^4 = 10000 and has a faithful representation of dimension 4,
# which can be extended to a faithful 10-dimensional representation.
# Therefore, the minimal kernel order of 1 is achievable.

# The largest possible dimension is the group order divided by the minimal kernel order.
max_dimension = group_order / min_kernel_order

# Print the final equation and the result.
print("The calculation for the largest possible dimension is as follows:")
print(f"{group_order} / {min_kernel_order} = {int(max_dimension)}")
