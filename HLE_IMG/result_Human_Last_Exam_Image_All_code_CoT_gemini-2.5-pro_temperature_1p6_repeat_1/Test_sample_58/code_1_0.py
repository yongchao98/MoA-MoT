import sympy

# Define the symbols
r, theta = sympy.symbols('r, θ')
h = 2 * r * theta * sympy.cot(theta)

# Create the equation object
equation = sympy.Eq(sympy.Symbol('h'), h)

# Print the final equation as a string
# We manually construct the string to ensure the numbers are explicitly shown
# as requested.
print(f"h = 2 * r * θ * cot(θ)")