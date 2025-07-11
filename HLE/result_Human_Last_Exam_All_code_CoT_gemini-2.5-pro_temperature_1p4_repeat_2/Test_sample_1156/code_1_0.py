import sympy

# Define the symbol
x = sympy.Symbol('x')

# Define the density function
# rho(x) = 1/(4*sqrt(x)) + 1/2
a, b, c, d = 1, 4, 1, 2
rho = a / (b * sympy.sqrt(x)) + sympy.Rational(c, d)

# Print the equation
equation = sympy.Eq(sympy.Symbol('rho(x)'), rho)
print("The normalised density of the invariant measure is given by the equation:")
# The problem asks to output each number in the final equation.
# The numbers are 1, 4, 1, 2. Let's make sure they are printed clearly.
print(f"rho(x) = {a}/({b}*sqrt(x)) + {c}/{d}")
