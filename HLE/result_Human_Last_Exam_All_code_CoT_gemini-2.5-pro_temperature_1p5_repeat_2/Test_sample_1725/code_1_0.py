# Part A: Dimension of pi's fibers
# The connection C is identified with the elastic distortion tensor in 3D space.
# This tensor is a general 3x3 matrix, as it represents the gradient of the 3D displacement vector field.
# The space of 3x3 matrices is a vector space of dimension 3 * 3 = 9.
dim_fibers = 3 * 3

# Part B: Number of coefficients specifying E
# The energy E is a quadratic form in the dislocation density tensor D, which is a general rank-2 tensor.
# E must be invariant under the cubic symmetry group O_h.
# The number of coefficients is the number of independent quadratic invariants of a rank-2 tensor under O_h.
# This is equal to the number of irreducible representations in the decomposition of the 9D space of rank-2 tensors.
# The decomposition under O_h is: A_1g (trace) + T_1g (antisymmetric part) + E_g + T_2g (symmetric traceless part).
# These are 4 distinct irreducible representations.
num_coefficients = 4

# Print the results in the required format.
# The first number is the answer to part A.
# The second number is the answer to part B.
print(f"{dim_fibers} {num_coefficients}")
