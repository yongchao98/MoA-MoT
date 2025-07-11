# The minimum curvature cost is the computational cost of the matrix inversion
# in the NGD update. We are looking for the lowest possible order of complexity.
# The number of samples is n, and the layer dimension is d, with n < d.

# Through a series of steps including exploiting the Kronecker product structure
# of the Fisher matrix and applying the Woodbury matrix identity, the most
# expensive part of the update, which is an inversion of a d^2 x d^2 matrix,
# can be reduced to an inversion of an n x n matrix.

# The cost of inverting an n x n matrix using standard algorithms like
# Gaussian elimination is O(n^3).

# Since n < d, this O(n^3) cost is the minimum achievable cost compared to
# other methods like naive inversion O(d^6) or a partial simplification O(d^3).

# The final equation for the order of complexity is n^3.
# The dominant part of this cost expression that determines its polynomial order
# is the exponent.

base_variable = 'n'
exponent_value = 3

print("The minimum achievable curvature cost has a polynomial order of complexity.")
print(f"The complexity is O({base_variable}^{exponent_value}).")
print("The final equation for the order of complexity is:")
print(f"Base: {base_variable}")
print(f"Exponent: {exponent_value}")