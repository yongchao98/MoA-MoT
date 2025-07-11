from sympy import symbols, integrate, pi, cos, sin

# Define the symbols for our coordinates
r, theta, z = symbols('r theta z')

# The function to integrate in cylindrical coordinates is f = z^2 * r^2.
# The volume element dV is r*dr*d(theta)*dz.
# So, the integrand is (z^2 * r^2) * r = z^2 * r^3.
integrand = z**2 * r**3

# Define the limits of integration for the cone.
# Radius r goes from 0 to the edge of the cone, which is r = 3/2 * (2-z).
# Height z goes from the base at 0 to the vertex at 2.
# Angle theta goes from 0 to 2*pi for the full cone.
r_limit_upper = (3/2) * (2 - z)
limits = (
    (r, 0, r_limit_upper),
    (z, 0, 2),
    (theta, 0, 2 * pi)
)

# Compute the triple integral
result = integrate(integrand, limits[0], limits[1], limits[2])

# The result is of the form (A * pi) / B
# We can extract the numerator and denominator to display them
# as requested.
p, q = result.as_numer_denom()
# The numerator contains pi, so we can divide by it to get the number
numerator = p / pi
denominator = q

print(f"The integral is calculated as: ({numerator} * pi) / {denominator}")
print("Final value of the integral:")
print(result)
