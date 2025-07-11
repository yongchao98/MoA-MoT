import sympy

# Define the symbols used in the equations
V0, r, pi, sigma1, sigma2 = sympy.symbols('V_0 r pi sigma_1 sigma_2')

# Numerators for the electric field expressions
# The common factor is 2 * V0
common_numerator_factor = 2 * V0

# Denominator is common for both E1 and E2
denominator = r * pi * (sigma1 + sigma2)

# Expression for E1 in region 1 (0 < phi < pi/2)
# The numerator includes sigma_2
numerator_E1 = common_numerator_factor * sigma2
E1 = numerator_E1 / denominator

# Expression for E2 in region 2 (pi/2 < phi < pi)
# The numerator includes sigma_1
numerator_E2 = common_numerator_factor * sigma1
E2 = numerator_E2 / denominator

# Print the derived expressions for the electric fields
# Using sympy.pretty_print for a more readable mathematical format
print("Electric field in Region 1 (E1):")
sympy.pretty_print(E1)
print("\nRepresents: E_1 = (2 * sigma_2 * V_0) / (r*pi*(sigma_1 + sigma_2)) * phi_hat")

print("\nElectric field in Region 2 (E2):")
sympy.pretty_print(E2)
print("\nRepresents: E_2 = (2 * sigma_1 * V_0) / (r*pi*(sigma_1 + sigma_2)) * phi_hat")