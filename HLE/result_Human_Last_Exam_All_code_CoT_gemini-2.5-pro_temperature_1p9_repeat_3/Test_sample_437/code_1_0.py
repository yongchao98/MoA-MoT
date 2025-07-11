import sympy as sp

# Define zeta as a positive real symbol
zeta = sp.symbols('zeta', real=True, positive=True)

# pi is a constant in sympy
pi = sp.pi

# Part 1: The square of the normalization constant for a 1s STO
# N = sqrt(zeta**3 / pi), so N^2 = zeta**3 / pi
N_squared = zeta**3 / pi

# Part 2: The integral over the angular coordinates (solid angle)
# Integral(d(phi) from 0 to 2*pi) * Integral(sin(theta)d(theta) from 0 to pi)
# = (2*pi) * (2) = 4*pi
angular_integral_result = 4 * pi

# Part 3: The integral over the radial part
# This is of the form integral(r^n * exp(-a*r))dr from 0 to infinity
# with n=1 and a=2*zeta. The standard result is n! / a^(n+1).
# So, the result is 1! / (2*zeta)^2 = 1 / (4*zeta**2)
radial_integral_result = 1 / (4 * zeta**2)

# Combine the parts to get the final result
final_result = N_squared * angular_integral_result * radial_integral_result

# --- Output the result ---
# The prompt requests that we output each number in the final equation.
# We will construct and print the equation step-by-step.
print("The integral is evaluated as the product of three parts:")
print(f"1. Squared Normalization Constant (N^2): {N_squared}")
print(f"2. Angular Integral (∫dΩ): {angular_integral_result}")
print(f"3. Radial Integral (∫ r*exp(-2ζr) dr): {radial_integral_result}")
print("\nThe final equation is:")
print(f"⟨ϕ_1s | 1/r | ϕ_1s⟩ = ({N_squared}) * ({angular_integral_result}) * ({radial_integral_result})")

# Print the simplified result
print("\nSimplifying the expression gives:")
# Use sp.simplify to ensure the expression is cancelled down correctly
simplified_result = sp.simplify(final_result)
print(f"⟨ϕ_1s | 1/r | ϕ_1s⟩ = {simplified_result}")