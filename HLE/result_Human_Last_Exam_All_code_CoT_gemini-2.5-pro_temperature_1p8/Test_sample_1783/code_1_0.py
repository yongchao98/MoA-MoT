import sympy
from sympy import pi, integrate, Rational

# 1. Define the symbolic variables
r, theta, z = sympy.symbols('r theta z')

# 2. Define the function and the integrand in cylindrical coordinates
# Original function f(x, y, z) = z^2 * (x^2 + y^2)
# In cylindrical coordinates, x^2 + y^2 = r^2, so f(r, theta, z) = z^2 * r^2
# The volume element dV is r*dr*d(theta)*dz
# The full expression to integrate is (z^2 * r^2) * r
integrand = z**2 * r**3

# 3. Define the limits of integration for the cone
# Height H = 2, Base Radius R = 3
# The radius of the cone at a given height z is r_cone = R * (1 - z/H)
r_limit_upper = 3 * (1 - z/Rational(2)) # Use Rational for precision
z_limit_lower = 0
z_limit_upper = 2
theta_limit_lower = 0
theta_limit_upper = 2 * pi

# 4. Calculate the integral
# We perform the integration from the inside out: dr, then dz, then d(theta)
# The result of the first two integrations (dr, dz) will be a constant value
# which we then multiply by the integral over theta (which is 2*pi).
integral_r = integrate(integrand, (r, 0, r_limit_upper))
integral_z = integrate(integral_r, (z, z_limit_lower, z_limit_upper))
final_integral = integrate(integral_z, (theta, theta_limit_lower, theta_limit_upper))

# 5. Display the results in a step-by-step manner
print("Problem: Find the integral of f(x,y,z) = z^2 * (x^2 + y^2) inside a cone of radius 3 and height 2.")
print("\nStep 1: The integral is set up in cylindrical coordinates.")
print(f"The integrand is f * r = (z^2 * r^2) * r = {integrand}")
print(f"The limits are:")
print(f"  - r from 0 to {r_limit_upper}")
print(f"  - z from {z_limit_lower} to {z_limit_upper}")
print(f"  - theta from {theta_limit_lower} to {theta_limit_upper}")

print("\nStep 2: The integral is solved.")
print(f"The integral of the expression with respect to r and z is:")
print(f"  Result(r, z) = {integral_z}")

print("\nStep 3: The final result is obtained by integrating with respect to theta.")
print("The final calculation is the result from Step 2 multiplied by the range of theta (2*pi).")
print("\n--- Final Equation ---")

# Extract numbers for the final equation printout
coeff = integral_z
angle = 2
angle_symbol = "pi"

# Display each number in the final equation
print(f"Value = ( {coeff} ) * ( {angle} * {angle_symbol} )")
print(f"Value = {final_integral}")

# Print the numerical approximation
print(f"\nNumerical value: {final_integral.evalf()}")

<<<108*pi/35>>>