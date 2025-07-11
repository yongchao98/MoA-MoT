import sympy as sp

# This script evaluates the integral <phi_1s| 1/r |phi_1s> for a 1s Slater-type orbital.
# The integral separates into three parts: a constant factor, an angular integral, and a radial integral.

# 1. Define the necessary symbols.
# r: radial distance
# theta: polar angle
# phi: azimuthal angle
# zeta: orbital exponent (a positive real number for a bound state)
r, theta, phi = sp.symbols('r theta phi')
zeta = sp.Symbol('zeta', positive=True)

print("Evaluating the integral <phi_1s| 1/r |phi_1s> = C * (Angular Integral) * (Radial Integral)")
print("--------------------------------------------------------------------------------------")

# 2. Define the constant part.
# The squared wavefunction is [ (zeta^3/pi)^(1/2) * exp(-zeta*r) ]^2 = (zeta^3 / pi) * exp(-2*zeta*r).
# The constant C is (zeta^3 / pi).
constants = zeta**3 / sp.pi
print(f"The constant C from the squared wavefunction is: {constants}")

# 3. Evaluate the angular part of the integral.
# The integral is over the solid angle, with the integrand sin(theta).
# Integral from phi=0 to 2*pi of Integral from theta=0 to pi of sin(theta) dtheta dphi.
angular_integrand = sp.sin(theta)
# Sympy performs the double integral.
angular_result = sp.integrate(angular_integrand, (phi, 0, 2 * sp.pi), (theta, 0, sp.pi))
print(f"The angular integral, Integral(sin(theta) dphi dtheta), evaluates to: {angular_result}")

# 4. Evaluate the radial part of the integral.
# The integrand from the wavefunction is exp(-2*zeta*r).
# The operator is 1/r.
# The volume element contributes r^2.
# So, the radial integrand is exp(-2*zeta*r) * (1/r) * r^2 = r * exp(-2*zeta*r).
radial_integrand = r * sp.exp(-2*zeta*r)
radial_result = sp.integrate(radial_integrand, (r, 0, sp.oo))
print(f"The radial integral, Integral(r * exp(-2*zeta*r) dr), evaluates to: {radial_result}")

# 5. Calculate the final result by multiplying the parts.
# Full Integral = C * (Angular Result) * (Radial Result)
final_result = constants * angular_result * radial_result

print("\nCombining these parts gives the final equation:")
# The following line shows how each calculated number/component contributes to the result.
print(f"<phi_1s| 1/r |phi_1s> = ({constants}) * ({angular_result}) * ({radial_result})")

# Print the simplified final answer.
final_result_simplified = sp.simplify(final_result)
print(f"\nAfter simplification, the final result is:")
print(f"<phi_1s| 1/r |phi_1s> = {final_result_simplified}")
