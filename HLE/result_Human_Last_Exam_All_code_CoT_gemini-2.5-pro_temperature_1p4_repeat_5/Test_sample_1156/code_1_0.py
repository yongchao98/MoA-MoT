import sympy

# Define the variable and the function
x = sympy.Symbol('x')
a = sympy.Rational(2, 3)
b = sympy.Rational(-1, 3)
h = a / sympy.sqrt(x) + b

# Print the final equation for the density
print("The normalised density of the invariant measure is h(x):")
# The user asked to output each number in the final equation.
# We will print the equation in a formatted way.
equation_str = f"h(x) = ({a}) / sqrt(x) + ({b})"
print(equation_str)
