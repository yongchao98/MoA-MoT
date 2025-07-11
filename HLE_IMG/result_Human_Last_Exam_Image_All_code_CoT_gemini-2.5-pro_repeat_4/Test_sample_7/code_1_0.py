import sympy

# Define the symbols for the formula
V2, mu1, mu2, e = sympy.symbols('V_2 μ_1 μ_2 e')

# The formula for the second plateau
# V_2 = (2 * mu_2 - mu_1) / e
# Coefficients
coeff_mu2 = 2
coeff_mu1 = -1

# Print the formula step-by-step as requested
print("The formula for the second voltage plateau (V_2) is:")
print(f"V_2 = ({coeff_mu2} * μ_2 - {abs(coeff_mu1)} * μ_1) / e")
