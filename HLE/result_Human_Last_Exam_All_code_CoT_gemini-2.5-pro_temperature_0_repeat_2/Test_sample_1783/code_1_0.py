import numpy as np
from scipy import integrate

# The function to integrate in cylindrical coordinates is f(r, theta, z) = z^2 * r^2.
# The volume element dV is r * dr * d(theta) * dz.
# So the integrand is (z^2 * r^2) * r = z^2 * r^3.
integrand = lambda r, theta, z: z**2 * r**3

# Define the limits of integration.
# z ranges from 0 to 2.
z_min, z_max = 0, 2

# theta ranges from 0 to 2*pi.
# These limits are functions of z.
theta_min = lambda z: 0
theta_max = lambda z: 2 * np.pi

# r ranges from 0 to (3/2)*(2-z).
# These limits are functions of z and theta.
r_min = lambda z, theta: 0
r_max = lambda z, theta: (3/2) * (2 - z)

# Use scipy's tplquad for triple integration.
# The order of integration is specified by the order of the limits: z, then theta, then r.
# tplquad(func, z_min, z_max, theta_min, theta_max, r_min, r_max)
result, error = integrate.tplquad(integrand, z_min, z_max, theta_min, theta_max, r_min, r_max)

# The analytical result is (108 * pi) / 35
numerator = 108
denominator = 35

print(f"The analytical result of the integral is the equation: ({numerator} * pi) / {denominator}")
print(f"The numerical value is: {result}")
print(f"The analytical value is: {numerator * np.pi / denominator}")
