import numpy as np
from scipy import integrate

# The function f(x,y,z) in cylindrical coordinates is z^2 * r^2.
# The volume element dV in cylindrical coordinates is r*dr*d(theta)*dz.
# So the integrand is (z^2 * r^2) * r = z^2 * r^3.
# The arguments to the lambda function must be in the order of integration:
# innermost to outermost. We integrate dr, then d(theta), then dz.
integrand = lambda r, theta, z: z**2 * r**3

# The integration ranges are defined from innermost to outermost.
# r_range: depends on z.
# theta_range: is constant.
# z_range: is constant.
# The ranges are specified for nquad as a list of functions of the outer variables.
# So r_range is a function of theta and z. theta_range is a function of z.
r_range = lambda theta, z: [0, 3 * (1 - z / 2)]
theta_range = lambda z: [0, 2 * np.pi]
z_range = [0, 2]

# Perform the numerical integration
result, error = integrate.nquad(integrand, [r_range, theta_range, z_range])

# The analytical result is (108 * pi) / 35
num_coeff = 108
den = 35

print("The final equation for the integral is (A * pi) / B")
print(f"A (Numerator coefficient): {num_coeff}")
print(f"B (Denominator): {den}")
print(f"Analytically, the integral is ({num_coeff} * pi) / {den}")
print(f"Numerical value: {result}")
print(f"Asserting analytical and numerical results are close: {np.isclose(result, (num_coeff * np.pi) / den)}")
