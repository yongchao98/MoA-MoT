import numpy as np
from scipy.integrate import quad

# Step 1: Define the orbital exponent zeta.
# For a hydrogen atom, zeta is 1.0. We will use this value for our calculation.
zeta = 1.0

# The integral we need to evaluate is of the form:
# I = 4 * zeta^3 * integral_from_0_to_inf(r * exp(-2 * zeta * r)) dr
# Let's define the function inside the integral (the integrand).
integrand = lambda r: r * np.exp(-2 * zeta * r)

# Step 2: Use scipy.integrate.quad to compute the radial integral.
# quad returns a tuple: (result, estimated_error)
radial_integral_value, error = quad(integrand, 0, np.inf)

# Step 3: Calculate the final result by plugging in all the numbers.
# The full expression is 4 * zeta^3 * (the radial integral we just computed).
prefactor = 4 * zeta**3
final_result = prefactor * radial_integral_value

# Step 4: Print out the components of the calculation and the final answer.
# This shows each number used in the final equation.
print(f"The calculation is for a 1s Slater orbital with exponent zeta = {zeta}")
print(f"The full integral to evaluate is: 4 * (zeta^3) * [Integral from 0 to infinity of r * exp(-2 * zeta * r) dr]")
print(f"The value of the radial integral is: {radial_integral_value:.4f}")
print("The final equation is: final_result = 4 * (zeta^3) * (radial_integral)")
print(f"Plugging in the numbers: {final_result:.4f} = {prefactor:.4f} * {radial_integral_value:.4f}")
print(f"\nThe final evaluated result is: {final_result:.4f}")
print(f"As expected from the analytical derivation, the result is equal to zeta ({zeta}).")

# Final answer block
# The question is symbolic, but the numerical evaluation for zeta=1.0 gives 1.0
final_answer_value = f"<<<{final_result:.1f}>>>"