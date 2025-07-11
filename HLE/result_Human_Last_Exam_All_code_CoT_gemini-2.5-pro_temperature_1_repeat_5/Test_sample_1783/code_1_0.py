import numpy as np
from scipy import integrate

# The problem is to find the integral of f(x,y,z) = z^2 * (x^2 + y^2)
# over a cone with base radius 3 and height 2.
# We convert to cylindrical coordinates (r, theta, z).
# The function becomes f = z^2 * r^2.
# The volume element dV is r * dr * d(theta) * dz.
# The integrand is thus (z^2 * r^2) * r = z^2 * r^3.

# We will use scipy.integrate.nquad, which integrates in the order of the provided
# function arguments. Our function will be func(r, z, theta), so the integration
# order will be dr, then dz, then d(theta).
integrand = lambda r, z, theta: z**2 * r**3

# Define the integration limits.
# The limit for r depends on z.
# The limits for z and theta are constant.
ranges = [
    [0, lambda z, theta: 3 * (1 - z/2)],  # r ranges from 0 to 3*(1-z/2)
    [0, 2],                               # z ranges from 0 to 2
    [0, 2 * np.pi]                        # theta ranges from 0 to 2*pi
]

# Perform the numerical integration
result, error = integrate.nquad(integrand, ranges)

# The analytical solution gives the exact result 108*pi/35.
# We will print the components of this exact result as requested.
numerator = 108
denominator = 35
analytical_value = numerator * np.pi / denominator

print(f"The integral of f(x,y,z)=z^2(x^2+y^2) over the given cone is calculated.")
print(f"The integral in cylindrical coordinates is:")
print(f"Integral from theta=0 to 2*pi, Integral from z=0 to 2, Integral from r=0 to 3(1-z/2) of (z^2 * r^3) dr dz dtheta")
print(f"The final equation for the exact result is ({numerator} * pi) / {denominator}")
print(f"Value from analytical solution: {analytical_value:.6f}")
print(f"Value from numerical integration: {result:.6f}")