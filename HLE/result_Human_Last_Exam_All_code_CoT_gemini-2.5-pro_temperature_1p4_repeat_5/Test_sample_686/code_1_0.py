import sympy

# Define symbolic variables for the parameters
beta, c, J = sympy.symbols('β c J')

# The expression for susceptibility derived from the steps above
# C_l = tanh(β*J)**l
# χ = β * Σ_{l=1 to ∞} c*(c-1)**(l-1) * tanh(β*J)**l
# χ = β * c * tanh(β*J) * Σ_{k=0 to ∞} [(c-1)*tanh(β*J)]**k
# The sum is a geometric series 1 / (1 - (c-1)*tanh(β*J))
numerator = beta * c * sympy.tanh(beta * J)
denominator = 1 - (c - 1) * sympy.tanh(beta * J)
chi_expression = numerator / denominator

# The note mentions N = β*c*(1-m_0**2)/(c-1).
# In the paramagnetic phase, m_0 = 0, so N = β*c/(c-1).
# We can express β*c as N*(c-1).
# chi_expression_with_N = N*(c-1)*sympy.tanh(beta*J) / (1-(c-1)*sympy.tanh(beta*J))

# The final output should contain each element of the equation.
# We will print the derived expression for chi.
# To make it clear, we can use a formatted string.

# Constructing the output string
# The instruction is to output each number, which we interpret as outputting the full formula.
final_equation_str = f"χ = ({beta} * {c} * tanh({beta}*{J})) / (1 - ({c} - 1) * tanh({beta}*{J}))"

print("The derived magnetic susceptibility in the paramagnetic phase is:")
print(final_equation_str)
