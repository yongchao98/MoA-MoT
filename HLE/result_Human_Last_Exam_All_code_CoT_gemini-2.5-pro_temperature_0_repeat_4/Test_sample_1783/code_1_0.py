import sympy
from sympy import pi, integrate

# Define the symbols for cylindrical coordinates
r, theta, z = sympy.symbols('r theta z')

# The function to integrate is f(x,y,z) = z^2 * (x^2 + y^2).
# In cylindrical coordinates, this becomes f = z^2 * r^2.
# The volume element dV is r * dr * d(theta) * dz.
# The full integrand is the function multiplied by the volume element factor 'r'.
integrand = z**2 * r**3

# The integration domain is a cone with base radius R=3 and height H=2.
# The base is on the z=0 plane, and the vertex is at (0,0,2).
# The limits of integration are:
# theta: 0 to 2*pi (a full circle)
# z: 0 to 2 (from base to vertex)
# r: 0 to 3*(1 - z/2). This is the radius of the cone at a given height z.

# We compute the definite triple integral.
# sympy.integrate can handle multiple integrals. The order of integration
# is from the inside out, so we list the variables in that order: r, then z, then theta.
result = integrate(integrand, (r, 0, 3 - 3*z/2), (z, 0, 2), (theta, 0, 2*pi))

# The final equation for the integral's value is 108 * pi / 35.
# Printing the symbolic result will display all the numbers in the final expression.
print("The exact value of the integral is:")
print(result)

# For clarity, we can also show the numerical approximation.
print("\nNumerical approximation:")
print(result.evalf())