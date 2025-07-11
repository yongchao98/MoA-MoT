import sympy

# Let p be the variable we are solving for.
p = sympy.Symbol('p')

# The two expressions for the exponent alpha as a function of p
# that are relevant for p > 4.
# alpha_2D comes from 2D Kakeya-type constructions.
# alpha_1D comes from 1D Kakeya-type constructions along the cone's rulings.
alpha_2D = 1/2 - 2/p
alpha_1D = 1/4 - 1/2/p

# The critical exponent occurs where the two behaviors are balanced,
# i.e., where the two alpha expressions are equal. This is where the
# slope of alpha(1/p) changes.
# We create an equation to find this critical point.
equation = sympy.Eq(alpha_2D, alpha_1D)

# Solve the equation for p.
# The result is a list of solutions. In this case, there is only one.
solution = sympy.solve(equation, p)
critical_exponent = solution[0]

# The equation to solve is:
# 1/2 - 2/p = 1/4 - 1/(2p)
# We can represent the numbers in the equation for clarity in the output.
a, b = 1/2, -2
c, d = 1/4, -1/2

print(f"The first expression for the exponent alpha is: {a} + ({b})/p")
print(f"The second expression for the exponent alpha is: {c} + ({d})/p")
print(f"We solve the equation: {a} + ({b})/p = {c} + ({d})/p")
print(f"The other critical exponent is p = {critical_exponent}")
