# Plan:
# 1. Determine the dimension of the fiber for the bundle pi.
#    The physical field C is the distortion tensor, a 3x3 matrix in 3D space.
#    The dimension of this space is calculated as dimension * dimension.
# 2. Determine the number of coefficients for the energy functional E.
#    E is a quadratic function of the strain tensor for a cubic crystal.
#    The number of independent elastic constants for a cubic crystal is a known value from solid-state physics.
# 3. Print the two results in the required format.

# Part A: Fiber Dimension
# The physical space is 3-dimensional.
spatial_dimension = 3
# The distortion tensor C is a rank-2 tensor in this space.
# The fiber of the bundle is the space of all possible values of C at a point,
# which is the space of 3x3 matrices.
fiber_dimension = spatial_dimension * spatial_dimension

# Part B: Number of Energy Coefficients
# The energy E for a cubic crystal is determined by its elastic constants.
# Due to the high symmetry of the cubic lattice, the 81 components of the
# general elasticity tensor reduce to just 3 independent coefficients.
# These are typically denoted c11, c12, and c44.
num_coefficients = 3

# Output the final answer in the format "dimension coefficients"
print(fiber_dimension, num_coefficients)
