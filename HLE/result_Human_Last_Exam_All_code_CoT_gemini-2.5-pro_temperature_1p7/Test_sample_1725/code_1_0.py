# Part A: Calculation of the dimension of the bundle's fibers.
# The connection C is identified with the distortion tensor in continuum mechanics,
# which is a second-rank tensor. At each point, it's a 3x3 matrix.
# The dimension of the space of 3x3 matrices is 3 * 3 = 9.
dimension_of_fiber = 9

# Part B: Calculation of the number of coefficients for the energy functional E.
# The energy E is a quadratic form in the dislocation density tensor alpha,
# where alpha is the curl of C.
# The number of coefficients is the number of independent quadratic invariants of a
# second-rank tensor under cubic symmetry.
# Using group representation theory, the space of second-rank tensors (dimension 9)
# decomposes into 4 distinct irreducible representations for the cubic group.
# Each irrep corresponds to one independent quadratic invariant.
# Therefore, there are 4 coefficients.
number_of_coefficients = 4

# The final answer is the dimension of the fiber followed by the number of coefficients.
# We print them in the requested format "A B".
print(f"{dimension_of_fiber} {number_of_coefficients}")