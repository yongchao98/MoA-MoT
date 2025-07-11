import sympy
import math

# Step 1: Define the symbols for our coordinates
r, theta, z = sympy.symbols('r theta z')

# Step 2: Define the integrand in cylindrical coordinates.
# The original function is f(x,y,z) = z^2 * (x^2 + y^2).
# In cylindrical coordinates, this is f(r,theta,z) = z^2 * r^2.
# The volume element dV is r*dr*dtheta*dz.
# The expression to integrate is f * r = z^2 * r^3.
integrand = z**2 * r**3

# Step 3: Define the limits of integration for the cone.
# The cone has height H=2 and base radius R=3.
# The radius of the cone at a height z is r(z) = R*(1 - z/H) = 3*(1 - z/2).
r_upper_limit = 3 * (1 - z/2)
z_lower_limit, z_upper_limit = 0, 2
theta_lower_limit, theta_upper_limit = 0, 2 * sympy.pi

# Step 4: Set up and print the integral expression.
print("The task is to find the integral of f(x,y,z) = z^2*(x^2+y^2) inside a cone of radius 3 and height 2.")
print("The integral in cylindrical coordinates is formulated as follows:")
print(f"  Integral from {theta_lower_limit} to {theta_upper_limit} d(theta) of")
print(f"    Integral from {z_lower_limit} to {z_upper_limit} d(z) of")
print(f"      Integral from 0 to {r_upper_limit} of ({integrand}) dr")
print("-" * 40)

# Step 5: Compute the integral using sympy.
integral_result = sympy.integrate(
    integrand,
    (r, 0, r_upper_limit),
    (z, z_lower_limit, z_upper_limit),
    (theta, theta_lower_limit, theta_upper_limit)
)

# Step 6: Display the final equation and its components.
# The result is of the form (A * pi) / B.
print(f"The exact result of the integration is: {integral_result}")

# Extract the numbers from the final equation as requested.
# The final equation is: Integral = 108 * pi / 35
# We can represent this as: Result = (NUMERATOR * pi) / DENOMINATOR
numerator_coeff, denominator = sympy.fraction(integral_result / sympy.pi)

print("\nIn the final equation, Result = (A * pi) / B:")
print(f"The number for the numerator coefficient, A, is: {numerator_coeff}")
print(f"The number for the denominator, B, is: {denominator}")

# Calculate the numerical value
numerical_value = float(integral_result)
print(f"\nThe approximate numerical value of the integral is: {numerical_value:.4f}")