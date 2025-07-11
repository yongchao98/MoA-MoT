import sympy

# The normalised invariant density is rho(x) = C * x^a
C = 0.5
a = -0.5

# We use sympy to pretty print the mathematical expression
x = sympy.Symbol('x')
rho = C * x**a

# Print the equation
print(f"The normalised density of the invariant measure is rho(x) = {sympy.pretty(rho)}")

# As requested, output each number in the final equation.
print(f"The equation for the density is rho(x) = C * x^a")
print(f"The coefficient C is: {C}")
print(f"The exponent a is: {a}")