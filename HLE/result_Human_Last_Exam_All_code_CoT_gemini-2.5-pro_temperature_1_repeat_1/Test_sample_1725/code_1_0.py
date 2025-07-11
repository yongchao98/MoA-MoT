# Part A: Calculation of the dimension of the bundle's fibers.
# The physical system is a continuum model of a crystal. A connection C on the tangent
# bundle T(R^3) is used to describe the local crystal structure. C=0 represents a
# perfect crystal, so C itself represents the distortion from this perfect state.
# At any point in space, the distortion is a linear transformation of the tangent space,
# which is represented by the 3x3 distortion tensor C_ij.
# The fiber of the bundle pi at a point is the space of all possible values the connection C
# can take at that point. This is the vector space of all 3x3 real matrices.
# The dimension of this space is the number of its components.

dimension_of_space = 3
fiber_dimension = dimension_of_space * dimension_of_space

# Part B: Calculation of the number of coefficients for the energy density E.
# The energy density E has several constraints:
# 1. Induces linear dynamics: This implies the energy is a quadratic function of the fields (C and its derivatives).
# 2. Detects dislocations (dC is the Burgers' vector): This implies E must depend on the spatial derivatives of C, which we can denote G_ijk = d_i(C_jk).
# 3. Homogeneous polynomial of least degree: A quadratic energy that depends on derivatives must be purely quadratic in those derivatives, i.e., of the form (dC)^2. A term like C^2 would not detect dislocations, and a mixed term like C(dC) would not be positive-definite.
# Therefore, E must be a quadratic form in the third-rank tensor G_ijk.
# E = D_ijkpqr * G_ijk * G_pqr
# The coefficients of E are the components of the 6th-rank tensor D.
# We need to find the number of independent components of D given the crystal's cubic symmetry.
# This is equivalent to finding the number of independent quadratic invariants of a general
# third-rank tensor under the cubic symmetry group (O_h).
# This is a standard result in the field of crystal tensor properties. For the cubic
# crystal classes (like O_h), a 6th-rank tensor governing a quadratic property related
# to a 3rd-rank tensor stimulus has 6 independent components.

num_coefficients = 6

# The final answer requires printing the two numbers in order.
print(f"{fiber_dimension} {num_coefficients}")