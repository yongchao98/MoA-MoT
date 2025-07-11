import sympy
from sympy import pi

# 1. Define the symbols for our coordinate system
r, theta, z = sympy.symbols('r theta z')

# 2. Define the integrand in cylindrical coordinates
# The function f(x,y,z) = z^2 * (x^2 + y^2) becomes f(r,theta,z) = z^2 * r^2.
# The volume element dV is r*dr*dθ*dz.
# The full integrand is f * r.
integrand = z**2 * r**3

# 3. Define the limits of integration for the cone
# The height 'z' goes from 0 to 2.
z_limits = (z, 0, 2)
# The angle 'theta' goes from 0 to 2*pi for the full cone.
theta_limits = (theta, 0, 2 * pi)
# The radius 'r' at a given height 'z' goes from 0 to the edge of the cone,
# which is described by the line r = 3 - (3/2)*z.
r_limits = (r, 0, 3 - (3/2)*z)

# 4. Compute the triple integral
# We integrate from the inside out: first 'r', then 'z', then 'theta'.
integral_value = sympy.integrate(integrand, r_limits, z_limits, theta_limits)

# 5. Output the final equation and its components
# The result is a fraction with pi. We extract the numerator and denominator.
if integral_value.has(pi):
    coeff = integral_value / pi
    numerator = sympy.numer(coeff)
    denominator = sympy.denom(coeff)

    print("The integral of f(x,y,z) = z^2 * (x^2 + y^2) over the specified cone is calculated.")
    print("\nThe final equation for the integral is:")
    print(f"Integral = ({numerator} * π) / {denominator}")
else:
    print(f"The result of the integral is: {integral_value}")
