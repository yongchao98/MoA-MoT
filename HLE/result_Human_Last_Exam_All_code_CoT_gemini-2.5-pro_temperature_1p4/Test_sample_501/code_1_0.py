import sympy

# Set up symbolic variables for the parameters in the problem.
# x: separation of the polymer ends
# l: length of a single strut (written as 'L' for clarity in sympy)
# n: number of mass points
# E0: kinetic energy of the polymer at zero extension
x, L, n, E0 = sympy.symbols('x l n E(0)')

# Define the derived force law expression.
# The formula is F = -(2 * E(0) * x / (n^2 * l^2)) * exp(x^2 / (n^2 * l^2))
Force = - (2 * E0 * x) / (n**2 * L**2) * sympy.exp(x**2 / (n**2 * L**2))

# Print the final force law in a readable format.
# The '2' is the number in the final equation.
print("The force law F(x) for the thermally isolated polymer is:")
sympy.pprint(Force, use_unicode=True)