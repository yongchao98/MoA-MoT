import sympy

# Define the symbols for the chemical potentials and charge
mu_2 = sympy.Symbol('μ_2')
mu_3 = sympy.Symbol('μ_3')
e = sympy.Symbol('e')
V_plateau = sympy.Symbol('V_plateau')

# The number involved in the equation
denominator_coefficient = 2

# Construct the formula
# The second voltage plateau corresponds to the transition between Stage 2 and Stage 3.
# The chemical potential of the plateau is approximated by the average of the
# chemical potentials of the two coexisting stages.
# Voltage is this average potential energy divided by the charge 'e'.
equation = (mu_2 + mu_3) / (denominator_coefficient * e)

# Print the final formula clearly, showing each component
print("The formula for the second voltage plateau is:")
print(f"{V_plateau} = ({mu_2} + {mu_3}) / ({denominator_coefficient} * {e})")