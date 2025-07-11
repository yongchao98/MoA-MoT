import sympy

# Define the symbols
MA, MB, MG, MH = sympy.symbols('MA MB MG MH')

# The problem is a theorem in geometry.
# The relationship between the lengths of the segments is given by:
# MG - MH = MB - MA
# This can also be written as MG + MA = MH + MB, which shows a nice symmetry.

# Let's represent this relationship as an equation
equation = sympy.Eq(MG - MH, MB - MA)

# We are asked to express MG - MH in terms of MA and MB.
# The expression is simply MB - MA.
# Let's print the derivation in the equation form as requested.
# To make it clear, we print the relationship as an equation string.
print(f"{equation.lhs} = {equation.rhs}")
