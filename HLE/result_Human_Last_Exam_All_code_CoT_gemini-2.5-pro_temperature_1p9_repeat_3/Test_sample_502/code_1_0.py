import math

# For the group G = (Z/10Z)^4 acting on C^10, the ring of invariants
# R^G is a polynomial ring.
# The action is chosen to be faithful, for instance, by letting the four
# generators of G act on the first four variables, and trivially on the rest.
# g_1 -> diag(exp(2*pi*i/10), 1, ..., 1)
# g_2 -> diag(1, exp(2*pi*i/10), 1, ..., 1)
# etc. for g_3, g_4.
# The invariant generators for this action are x_1^10, x_2^10, x_3^10, x_4^10,
# and x_5, x_6, x_7, x_8, x_9, x_10.

# Degrees of the fundamental invariants
degrees = [10, 10, 10, 10, 1, 1, 1, 1, 1, 1]

# According to the Chevalley-Shephard-Todd theorem, for a reflection group,
# the dimension of the coinvariant algebra is the product of the degrees
# of the fundamental invariants. This dimension is also equal to the order of the group.
dimension = math.prod(degrees)

# Create the equation string
equation_parts = [str(d) for d in degrees]
equation_str = " * ".join(equation_parts)

print("The degrees of the fundamental invariants are:", degrees)
print(f"The dimension of R/I is the product of these degrees:")
print(f"{equation_str} = {dimension}")
print(f"\nThus, the largest possible dimension for R/I is {dimension}.")
