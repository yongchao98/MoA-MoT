# Part A: Dimension of the fibers of the bundle pi
# The connection C represents the distortion field, which is a map from the tangent space T(R^3)
# to an internal vector space (also R^3). This map is represented by a 3x3 matrix at each point.
# The number of components in this matrix is the dimension of the fiber.
dim_fiber = 3 * 3

# Part B: Number of coefficients specifying the energy E
# The energy E is a quadratic function of the dislocation density tensor alpha (dC).
# We need to find the number of independent quadratic invariants of a 3x3 tensor alpha
# under cubic crystal symmetry.
# We decompose alpha into its symmetric part (S) and antisymmetric part (A).
# Number of invariants for the symmetric part S (classical elasticity for cubic crystals) = 3.
# Number of invariants for the antisymmetric part A (a pseudovector) = 1.
# The total number of coefficients is the sum.
num_coeffs = 3 + 1

# Output the final answer in the specified format "A B"
print(f"{dim_fiber} {num_coeffs}")