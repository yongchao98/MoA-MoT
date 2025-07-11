import sympy as sp

# Define the variables for cylindrical coordinates
r, theta, z = sp.symbols('r theta z')

# Define the function f(x,y,z) = z^2 * (x^2 + y^2) in cylindrical coordinates
# x^2 + y^2 = r^2, so f = z^2 * r^2
f_cylindrical = z**2 * r**2

# The integrand includes the Jacobian 'r' for the volume element dV = r*dr*d_theta*dz
integrand = f_cylindrical * r

# Define the limits of integration for the cone
# Height H = 2, Radius R = 3
# Limit for theta: 0 to 2*pi
theta_limits = (theta, 0, 2*sp.pi)
# Limit for z: 0 to H=2
z_limits = (z, 0, 2)
# Limit for r: 0 to R*(1 - z/H) = 3*(1 - z/2)
r_upper_limit = 3 * (1 - z/2)
r_limits = (r, 0, r_upper_limit)

# Calculate the triple integral
integral_result = sp.integrate(integrand, r_limits, z_limits, theta_limits)

# --- Output the results ---
print("To find the integral of f(x,y,z) = z^2*(x^2+y^2) over the specified cone, we solve the following triple integral in cylindrical coordinates:")
print("\nIntegral Equation:")
print(f"  I = Integral from 0 to 2*pi d(theta) * Integral from 0 to 2 d(z) * Integral from 0 to 3*(1-z/2) d(r) of ({integrand})")

# Deconstruct the result for clear printing
# The result from sympy is of the form (numerator/denominator) * pi
# E.g., 108*pi/35
num, den = integral_result.as_numer_denom()
# The numerator will contain pi, so we can separate it
coeff = num / sp.pi

print("\nFinal Equation with Result:")
print(f"  I = {int(coeff)} * pi / {int(den)}")
print(f"\nThe numerical value is approximately {integral_result.evalf()}")
