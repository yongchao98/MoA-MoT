# Number of variables in the integral I(a)
D_ambient = 9

# Dimension of the integration domain (dx dy)
D_integration = 2

# The degree of the polynomial phase P(x, y; a)
deg_poly = 3

# To make the phase maximally flat, we set its Taylor coefficients for degrees less than deg_poly to zero.
# The number of coefficients for a 2D polynomial of degree k is (k+1)*(k+2)/2.
# We impose conditions on terms of degree 0, 1, and 2.
# We don't need to set the constant term (degree 0) to zero, as it only adds a phase factor.
num_coeffs_deg0 = 1
num_coeffs_deg1 = 2
num_coeffs_deg2 = 3

# Number of conditions to define the subspace of maximal flatness
num_conditions = num_coeffs_deg1 + num_coeffs_deg2
print(f"Number of conditions to achieve maximal flatness: {num_conditions}")

# Dimension of the subspace where the decay of I(a) is slowest
d_subspace = D_ambient - num_conditions
print(f"Dimension of the subspace of slowest decay: {d_subspace}")

# On this subspace, the phase is a homogeneous polynomial of degree kappa
kappa = deg_poly
print(f"Order of the phase function on this subspace: {kappa}")

# The decay exponent sigma for a D_integration-dimensional integral with a phase of order kappa
sigma = D_integration / kappa
print(f"Decay exponent sigma = D_integration / kappa = {D_integration}/{kappa} = {sigma}")

# The integral over the subspace diverges if p <= d_subspace / sigma
p_critical = d_subspace / sigma
print(f"The integral diverges for p <= d_subspace / sigma")
print(f"p <= {d_subspace} / {sigma}")
print(f"p <= {p_critical}")

print(f"\nThe largest value of p such that I is not in L^p(R^9) is {p_critical}.")
