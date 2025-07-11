import sympy as sp

# Plan:
# 1. State the problem: We are evaluating the expectation value <φ_1s| 1/r |φ_1s>.
# 2. Define the integral in mathematical terms using spherical coordinates.
# 3. Define the normalized 1s Slater-type orbital (STO) φ_1s.
# 4. Set up the full integral expression.
# 5. Separate the integral into radial and angular parts.
# 6. Use sympy to evaluate the angular part.
# 7. Use sympy to evaluate the radial part.
# 8. Combine the results to get the final answer and display the calculation.

print("This script evaluates the integral <φ_1s| 1/r |φ_1s> for a 1s Slater-type orbital (STO).")
print("This is the expectation value of the potential energy operator 1/r.")
print("-" * 50)

# --- Symbolic Setup ---
# Define symbols for the integration variables and the orbital exponent zeta
r, theta, phi = sp.symbols('r theta phi', real=True, positive=True)
zeta = sp.symbols('zeta', real=True, positive=True)

# The normalized 1s STO wavefunction is φ_1s = sqrt(zeta^3 / π) * exp(-zeta*r)
phi_1s = sp.sqrt(zeta**3 / sp.pi) * sp.exp(-zeta * r)

# The integral is I = ∫∫∫ (φ_1s*) * (1/r) * (φ_1s) * r^2 * sin(θ) dr dθ dφ
# Since φ_1s is real, φ_1s* = φ_1s.
integrand = phi_1s * (1/r) * phi_1s * r**2 * sp.sin(theta)

# --- Explanation and Separation ---
print(f"The normalized 1s Slater orbital is defined as: φ_1s(r) = {phi_1s}")
print("The operator is 1/r.")
print("The volume element is dτ = r^2 * sin(θ) dr dθ dφ.")
print("\nThe full integral is I = ∫∫∫ φ_1s * (1/r) * φ_1s * dτ over all space.")
print("This can be separated into a radial part and an angular part.")
print("I = ( ∫[0, ∞] [radial part] dr ) * ( ∫[0, 2π]dφ ∫[0, π] [angular part] dθ )")
print("-" * 50)


# --- Step-by-step Evaluation ---
# 1. Angular Part
angular_integral = sp.integrate(sp.sin(theta), (theta, 0, sp.pi)) * sp.integrate(1, (phi, 0, 2*sp.pi))

print("1. Evaluating the Angular Part:")
print(f"   ∫[0, 2π]dφ ∫[0, π]sin(θ)dθ = ({sp.integrate(1, (phi, 0, 2*sp.pi))}) * ({sp.integrate(sp.sin(theta), (theta, 0, sp.pi))})")
print(f"   Result of angular integral = {angular_integral}")
print("-" * 50)

# 2. Radial Part
radial_integrand = phi_1s * (1/r) * phi_1s * r**2
radial_integrand_simplified = sp.simplify(radial_integrand)
radial_integral = sp.integrate(radial_integrand_simplified, (r, 0, sp.oo))

print("2. Evaluating the Radial Part:")
print(f"   The radial integrand is: {radial_integrand_simplified}")
print(f"   ∫[0, ∞] ({radial_integrand_simplified}) dr = {radial_integral}")
print("-" * 50)

# 3. Final Combination
final_result = radial_integral * angular_integral

print("3. Combining the parts to get the final answer:")
print("   Final Value = (Result of radial part) * (Result of angular part)")
# The request is to "output each number in the final equation"
print(f"   Final Value = ({radial_integral}) * ({angular_integral}) = {sp.simplify(final_result)}")
