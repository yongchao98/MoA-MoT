import math

# This script derives and prints the analytical expression for the initial
# magnetization curve M(H) of a thin superconducting bar in a perpendicular
# magnetic field, based on the critical-state model.

# -----------
# Definitions
# -----------
# M: Magnetization along the y-axis (in A/m)
# H: Applied uniform magnetic field along the y-axis (in A/m)
# Jc: Critical current density, a constant (in A/m^2)
# a: Half-width of the bar's cross-section (along x-axis, in m)
# b: Half-thickness of the bar's cross-section (along y-axis, in m)
# The problem assumes a thin bar limit, where b << a.
# pi: Mathematical constant pi

# -----------
# The Final Derived Expression
# -----------
# The magnetization M is found to be:
# M(H) = - (Jc * a / 4) * tanh^2( (pi * H) / (2 * b * Jc) )
# where tanh is the hyperbolic tangent function.

# Let's break down the final equation and print it for clarity.
print("The analytical expression for the initial magnetization M as a function of the applied field H is composed of several parts:")
print("-" * 80)

# Define the components of the equation as string variables
term_M = "M(H)"
sign = "-"
coefficient_part = "(Jc * a / 4)"
function_part = "tanh^2" # hyperbolic tangent squared
argument_numerator = "(pi * H)"
argument_denominator = "(2 * b * Jc)"

# Print the full equation in a structured way
print(f"The equation is: {term_M} = {sign} {coefficient_part} * {function_part}[ {argument_numerator} / {argument_denominator} ]")
print("-" * 80)

# Print each part of the equation as requested
print("Component-wise breakdown:")
print(f"1. The magnetization: {term_M}")
print(f"2. The sign, indicating diamagnetism: '{sign}'")
print(f"3. The saturation magnetization coefficient: {coefficient_part}")
print(f"4. The function describing the field dependence: '{function_part}' (hyperbolic tangent squared)")
print(f"5. The argument of the function, which is a dimensionless ratio of fields:")
print(f"   - Numerator: {argument_numerator}")
print(f"   - Denominator: {argument_denominator}")
