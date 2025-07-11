import math

# This script evaluates the electron-nucleus attraction integral <phi|1/r|phi>
# for a 1s Slater-type orbital (STO). This integral represents the average
# potential energy of an electron in this orbital due to the nucleus.

# For this demonstration, we'll use a common orbital exponent for Helium.
zeta = 1.24

print("Evaluation of <phi_1s | 1/r | phi_1s>")
print("-" * 40)
print(f"The 1s Slater orbital is: phi(r) = N * exp(-zeta*r)")
print(f"The chosen orbital exponent is zeta = {zeta:.4f}\n")

print("The integral to evaluate is I = integral( (phi(r))^2 / r dV )\n")

# Step 1: Write the integral using the normalization constant N^2 = zeta^3 / pi
# and the volume element in spherical coordinates dV = r^2 sin(theta) dr d(theta) d(phi).
print("Step 1: The integral in spherical coordinates")
print("I = (zeta^3 / pi) * integral_0^inf( r * exp(-2*zeta*r) dr ) * integral_angles( d(Omega) )")
N_squared = zeta**3 / math.pi
print(f"The term N^2 = zeta^3 / pi = ({zeta**3:.4f}) / pi = {N_squared:.4f}\n")

# Step 2: Solve the angular part of the integral.
angular_integral_val = 4 * math.pi
print("Step 2: The angular integral")
print(f"The integral over solid angle d(Omega) is 4*pi, approx {angular_integral_val:.4f}.\n")

# Step 3: Solve the radial part of the integral.
print("Step 3: The radial integral")
print("The radial part is integral_0^inf( r * exp(-2*zeta*r) dr ).")
# Using the formula integral(x*e^(-ax))dx = 1/a^2, with a = 2*zeta
a = 2 * zeta
radial_integral_val = 1 / (a**2)
print(f"This evaluates to 1 / (2*zeta)^2 = 1 / ({a:.2f})^2 = {radial_integral_val:.4f}\n")

# Step 4: Combine all parts to get the final result.
print("Step 4: Combine all parts")
print("I = (N^2) * (Radial Part) * (Angular Part)")
print(f"I = ({N_squared:.4f}) * ({radial_integral_val:.4f}) * ({angular_integral_val:.4f})")
# The analytical result simplifies to just zeta.
final_value = zeta
print("This simplifies analytically to: I = zeta\n")

# Final result presentation
print("Final Result:")
print("The final equation for the integral is:")
print("I = zeta")
print(f"\nFor our chosen value of zeta = {zeta:.4f}, the integral evaluates to:")
print(f"I = {final_value:.4f}")
