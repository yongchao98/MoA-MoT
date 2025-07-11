import sympy

# Define the symbols
V_2 = sympy.Symbol('V_2')  # Voltage of the second plateau
mu_2 = sympy.Symbol('μ_2')  # Chemical potential for the second plateau transition
mu_Li = sympy.Symbol('μ_Li') # Chemical potential of the Li/Li+ reference electrode
e = sympy.Symbol('e')        # Elementary charge

# Create the equation
# The voltage V_2 is the negative of the difference between the chemical potential
# of Li in the graphite at stage 2 (μ_2) and the chemical potential of the
# Li reference electrode (μ_Li), all divided by the elementary charge (e).
equation = sympy.Eq(V_2, -(mu_2 - mu_Li) / e)

# Print the formula part of the equation
formula_rhs = equation.rhs

# To match the requested format, we print the formula side of the equation.
# The `sympy.printing.pretty` function makes it look nicer.
# We are asked to output the formula for the plateau.
print("The formula for the second plateau voltage (V_2) is:")
print(sympy.printing.pretty(formula_rhs, use_unicode=True))

# The problem asks to also output numbers, but the formula is symbolic.
# We can show how the expression relates to the given values.
# V_2 = 0.13 V
# V_1 = 0.09 V
# The difference in voltage between the plateaus can be expressed as:
# V_2 - V_1 = (μ_1 - μ_2) / e
# 0.13 - 0.09 = 0.04 V, so (μ_1 - μ_2) / e = 0.04 V.
# This does not give a formula for V_2 in terms of only μ's, but it's a key relationship.
# The fundamental formula remains the most direct answer.