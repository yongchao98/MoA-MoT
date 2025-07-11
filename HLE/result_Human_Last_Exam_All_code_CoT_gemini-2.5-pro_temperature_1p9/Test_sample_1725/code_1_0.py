# Plan:
# Part A: Determine the dimension of the fiber of the bundle pi.
# The connection C is identified with the crystal's distortion tensor. In 3D space, this
# is represented by a 3x3 matrix at each point. The fiber of the bundle pi is the space
# of these matrices. The dimension is the number of elements in a 3x3 matrix.

# Part B: Determine the number of coefficients for the energy density E.
# The energy E is a quadratic function of the strain tensor (the symmetric part of C).
# The coefficients form the stiffness tensor. For a "cubic, axis-aligned crystal",
# material symmetry constrains the number of independent elastic constants.

# Calculation for Part A
# The space is the set of all 3x3 real matrices.
dimension_of_space = 3
dimension_of_fiber = dimension_of_space * dimension_of_space

# Calculation for Part B
# For a cubic crystal system, the elastic stiffness tensor, which defines the energy
# functional, has 3 independent components (C11, C12, C44 in Voigt notation).
num_coefficients = 3

# Output the results in the required format.
print(f"{dimension_of_fiber} {num_coefficients}")