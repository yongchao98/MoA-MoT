import sympy

# Define the symbol x
x = sympy.Symbol('x')

# Define the symbolic constants from the graph
a, b, c, d = sympy.symbols('a b c d')

# Construct the numerator based on the x-intercepts
# Root at -b and b -> (x+b)(x-b) = x**2 - b**2
# Double root at d -> (x-d)**2
numerator = (x**2 - b**2) * (x - d)**2

# Construct the denominator based on the vertical asymptotes
# Asymptote at a (odd) -> (x-a)
# Asymptote at c (even) -> (x-c)**2
denominator = (x - a) * (x - c)**2

# Create the expression for the function f(x)
fx = numerator / denominator

# Create a nicely formatted string for the equation
# Using sympy.pretty_print would be an option, but for this specific format,
# manual string creation is more direct.
equation_str = f"f(x) = {sympy.printing.pretty(numerator, use_unicode=False)} / {sympy.printing.pretty(denominator, use_unicode=False)}"

# Manually replace the structure for a single line print
equation_str_oneline = f"f(x) = ((x**2 - b**2)*(x - d)**2) / ((x - a)*(x - c)**2)"

# Let's create a final string that shows all numbers like 2 in the powers clearly.
# Remember in the final code you still need to output each number in the final equation!
final_equation_string = f"f(x) = ((x^2 - b^2)(x - d)^2) / ((x - a)(x - c)^2)"

print("Based on the analysis of the graph's features, the equation for f(x) is:")
print(final_equation_string)
