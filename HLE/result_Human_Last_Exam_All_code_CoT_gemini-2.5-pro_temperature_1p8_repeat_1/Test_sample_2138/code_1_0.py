import math
import cmath

# The integral evaluates to a complex number.
# The analytical solution is I = i * (2*pi)/(3*sqrt(3)) * (Gamma(1/3)/Gamma(2/3))^3

# Calculate the components of the formula
gamma_1_3 = math.gamma(1/3)
gamma_2_3 = math.gamma(2/3)
pi_val = math.pi
sqrt3 = math.sqrt(3)

# Calculate the main coefficients
coeff1 = (2 * pi_val) / (3 * sqrt3)
ratio_gamma_cubed = (gamma_1_3 / gamma_2_3)**3

# Calculate the final value
# The result is purely imaginary, so we represent it as a complex number.
result = 1j * coeff1 * ratio_gamma_cubed

# Print the equation with its components
print(f"The integral is evaluated using the formula:")
print(f"I = i * (2*pi / (3*sqrt(3))) * (Gamma(1/3) / Gamma(2/3))^3")
print("\nWhere:")
print(f"  Gamma(1/3)          = {gamma_1_3}")
print(f"  Gamma(2/3)          = {gamma_2_3}")
print(f"  2*pi / (3*sqrt(3))    = {coeff1}")
print(f"  (Gamma(1/3)/Gamma(2/3))^3 = {ratio_gamma_cubed}")
print("\nFinal equation with intermediate values:")
# Note: "j" is Python's notation for the imaginary unit i.
print(f"I = {1j} * {coeff1} * {ratio_gamma_cubed}")

print("\nResult:")
# Using print() on a complex number will format it correctly.
print(result)
