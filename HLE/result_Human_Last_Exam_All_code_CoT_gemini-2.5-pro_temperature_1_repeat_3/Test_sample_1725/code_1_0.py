# Part A: Dimension of the fibers
# The physical field is the strain tensor, which is a symmetric 3x3 matrix.
# The dimension of the space of symmetric 3x3 matrices is (3 * (3 + 1)) / 2.
dim_A = (3 * 4) // 2

# Part B: Number of coefficients specifying E
# The energy E is a quadratic form in the derivative of the strain tensor, ε'.
# E must be invariant under the cubic crystal symmetry group (O_h).
# The number of independent coefficients in the quadratic form is equal to
# the number of independent elastic constants for a cubic crystal.
# This is a standard result from solid-state physics and elasticity theory.
num_coeffs_B = 3

# The problem asks for the answer in the format "A B".
# The final equation is simply the pair of results.
print(f"{dim_A} {num_coeffs_B}")