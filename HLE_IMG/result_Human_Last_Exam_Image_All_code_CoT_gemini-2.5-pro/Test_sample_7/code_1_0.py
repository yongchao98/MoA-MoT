import sympy

# Define the symbols for the formula
# V represents the voltage of the plateau
# mu represents the chemical potential
# k is the stage number
# e represents the elementary charge

# For the second plateau, the transition leads to the formation of Stage 2.
# Therefore, the relevant stage number is 2.
k = 2

# The voltage (V) is related to the chemical potential (μ_k) and the elementary charge (e)
# by the formula V = -μ_k / e.
V, mu, e = sympy.symbols('V μ e')
formula_eq = sympy.Eq(V, -sympy.Symbol(f'μ_{k}') / e)

# The question asks for the formula for the voltage, which is the right side of the equation.
formula_str = str(formula_eq.rhs)

# The prompt requires outputting each number in the final equation.
# In our formula "-μ_2 / e", the number is 2. The code explicitly uses this number.
print("A simple formula that best approximates the second plateau is:")
print(formula_str)
