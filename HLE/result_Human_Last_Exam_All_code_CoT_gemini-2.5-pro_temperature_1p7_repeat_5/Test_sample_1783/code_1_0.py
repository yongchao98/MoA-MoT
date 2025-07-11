import sympy

# Define the variables for integration
r, theta, z = sympy.symbols('r theta z')

# Define the function f in cylindrical coordinates
# Original function: f(x, y, z) = z^2 * (x^2 + y^2)
# In cylindrical coordinates: f(r, z) = z^2 * r^2
f_cylindrical = z**2 * r**2

# The volume element in cylindrical coordinates is r * dr * d(theta) * dz
# The complete integrand is f_cylindrical times r
integrand = f_cylindrical * r

# Define the limits of integration for the cone
# Height (H) = 2, Radius (R) = 3
# The radius 'r' at a given height 'z' is described by the line from (r=3, z=0) to (r=0, z=2).
# Equation: r = R * (1 - z/H) = 3 * (1 - z/2)
r_upper_limit = 3 * (1 - z/2)

# Set up and compute the triple integral.
# Order of integration is dr, then dz, then dtheta.
integral_result = sympy.integrate(
    integrand,
    (r, 0, r_upper_limit),
    (z, 0, 2),
    (theta, 0, 2 * sympy.pi)
)

# Print the final symbolic result. The components of this expression
# are the numbers in the final equation for the integral's value.
print("The exact value of the integral is:")
# The output 108*pi/35 contains the numbers 108, pi, and 35.
print(integral_result)