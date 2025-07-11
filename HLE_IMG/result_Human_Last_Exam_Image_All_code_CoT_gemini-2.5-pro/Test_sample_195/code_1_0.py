import sympy

# Define the symbolic variables
x, a, b, c, d = sympy.symbols('x a b c d')

# Construct the numerator based on the roots: -b, b, d
# The factors are (x - (-b)), (x - b), and (x - d)
numerator = (x**2 - b**2) * (x - d)

# Construct the denominator based on the vertical asymptotes: a, c
# The factors are (x - a) and (x - c)
denominator = (x - a) * (x - c)

# Create the function f(x)
# f_x = numerator / denominator
# We will print it in a more readable format.

# Output the final equation
# The problem asks to output each number in the equation.
# Since we have symbolic constants, we will format the string to show them clearly.
print(f"The equation for the function is:")
print(f"f(x) = ((x^2 - b^2) * (x - d)) / ((x - a) * (x - c))")
print("\nExpanded form:")
print(f"Numerator: {sympy.expand(numerator)}")
print(f"Denominator: {sympy.expand(denominator)}")
