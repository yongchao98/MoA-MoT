import sympy
from sympy import symbols, integrate, exp, pi, oo, sin, pretty

# Define symbolic variables
r = symbols('r', real=True, positive=True)
zeta = symbols('zeta', real=True, positive=True)
theta, phi = symbols('theta phi', real=True)

# --- Introduction ---
print("This script evaluates the integral: <phi_1s | 1/r | phi_1s>\n")
print("This integral represents the expectation value of the electron-nucleus attraction potential")
print("for an electron described by a 1s Slater-type orbital (STO).\n")

# --- Components of the Integral ---
print("--- 1. Components of the Integral ---")

N_squared_expr = zeta**3 / pi

print(f"The square of the normalized 1s STO is: phi_1s^2 = ({pretty(N_squared_expr, use_unicode=False)}) * exp(-2*zeta*r)")
print(f"The operator is: 1/r")
print(f"The volume element in spherical coordinates is: d_tau = r^2 * sin(theta) * dr * d(theta) * d(phi)\n")

# --- Setting up the Integral ---
print("--- 2. Setting up the Full Integral ---")
print("The integral I is the product of these terms integrated over all space:")
print("I = Integral( (zeta^3/pi) * exp(-2*zeta*r) * (1/r) * r^2 * sin(theta) ) dr d(theta) d(phi)")
print("This simplifies by combining r terms (r^2/r = r):")
print("I = (zeta^3/pi) * Integral( r * exp(-2*zeta*r) * sin(theta) ) dr d(theta) d(phi)\n")

# --- Solving the Integral ---
print("--- 3. Solving the Integral by Parts ---")
print("The integral can be separated into a constant factor, an angular part, and a radial part.\n")

# Angular Part
print("  a) Angular Integration:")
angular_integrand = sin(theta)
angular_integral_result = integrate(angular_integrand, (phi, 0, 2*pi), (theta, 0, pi))
print(f"     Integral( sin(theta) d(phi) d(theta) ) from phi=0..2*pi, theta=0..pi  =  {angular_integral_result}")

# Radial Part
print("\n  b) Radial Integration:")
radial_integrand_expr = r * exp(-2 * zeta * r)
radial_integral_result = integrate(radial_integrand_expr, (r, 0, oo))
radial_pretty = pretty(radial_integral_result, use_unicode=False)
print(f"     Integral( {pretty(radial_integrand_expr, use_unicode=False)} ) dr from r=0..oo  =  {radial_pretty}")

# --- Final Calculation ---
print("\n--- 4. Final Calculation ---")
print("Now, we combine all the parts: I = (Constant Factor) * (Angular Result) * (Radial Result)")

print("\nHere is the final equation with all the symbolic parts and numbers:")
# Using string formatting to show the equation as requested
# The numbers involved are: 3, 4, 1, 4, 2
term1_str = "zeta**3 / pi"
term2_str = "4*pi" 
term3_str = "1 / (4*zeta**2)"

print(f"I = ({term1_str}) * ({term2_str}) * ({term3_str})")

print("\nSimplifying the expression:")
print("- The 'pi' in the numerator and denominator cancel out.")
print("- The '4' in the numerator and denominator cancel out.")
print("- The remaining term is zeta**3 / zeta**2, which simplifies to zeta.")

final_result = N_squared_expr * angular_integral_result * radial_integral_result

print("\n" + "="*40)
print(f"The final simplified result of the integral is: {final_result}")
print("="*40)