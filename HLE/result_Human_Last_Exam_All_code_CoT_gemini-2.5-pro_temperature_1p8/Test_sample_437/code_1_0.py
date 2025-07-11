import sympy as sp

# Set up the symbols for our calculation
# r, theta, phi are the spherical coordinates
# zeta is the orbital exponent of the Slater orbital
r, theta, phi = sp.symbols('r theta phi', real=True, positive=True)
zeta = sp.symbols('zeta', real=True, positive=True)

# --- Step 1: Define the 1s Slater Orbital ---
# The normalization constant N for a 1s Slater orbital is (zeta^3 / pi)^(1/2)
N = sp.sqrt(zeta**3 / sp.pi)

# The 1s Slater orbital phi_1s = N * exp(-zeta * r)
phi_1s = N * sp.exp(-zeta * r)

print("--- The Problem ---")
print(f"We are evaluating the integral <phi_1s | 1/r | phi_1s> for the orbital:")
print(f"phi_1s = {phi_1s}\n")


# --- Step 2: Set up the integral ---
# The integral is Integral(phi_1s^* * (1/r) * phi_1s * d_tau)
# Since phi_1s is real, phi_1s^* = phi_1s.
# The volume element d_tau in spherical coordinates is r^2 * sin(theta) dr d(theta) d(phi)
# The operator is 1/r.

# Full integrand:
integrand = phi_1s * (1/r) * phi_1s * r**2 * sp.sin(theta)

print("--- The Setup ---")
print("The integrand in spherical coordinates is:")
# Using sp.pretty_print for a clearer mathematical output
sp.pretty_print(integrand)
print()


# --- Step 3: Solve the integral by separating it ---
# The integral can be separated into angular and radial parts:
# Integral = (N^2) * Integral_r(r * exp(-2*zeta*r) dr) * Integral_theta(sin(theta) d(theta)) * Integral_phi(d(phi))

# Calculate the angular parts:
integral_over_phi = sp.integrate(1, (phi, 0, 2*sp.pi))
integral_over_theta = sp.integrate(sp.sin(theta), (theta, 0, sp.pi))
angular_integral_result = integral_over_phi * integral_over_theta

print("--- Calculation Steps ---")
print(f"1. The normalization constant squared (N^2) is: {N**2}")
print(f"2. The angular integral part is: {angular_integral_result}")

# Calculate the radial part:
radial_integrand = r * sp.exp(-2*zeta*r)
radial_integral_result = sp.integrate(radial_integrand, (r, 0, sp.oo))

print(f"3. The radial integral part is: {radial_integral_result}\n")

# --- Step 4: Combine all parts to get the final result ---
# The total integral is the product of these three parts.
print("--- Final Equation ---")
print("The final value is the product of the three parts calculated above:")

# Here we output the components of the final equation as requested
equation_str = f"Value = ({N**2}) * ({angular_integral_result}) * ({radial_integral_result})"
print(equation_str)

final_result = sp.simplify(N**2 * angular_integral_result * radial_integral_result)

print("\n--- Final Result ---")
print(f"The simplified result of the integral <phi_1s | 1/r | phi_1s> is: {final_result}")
