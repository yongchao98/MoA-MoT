# This script formulates and prints the equation for the function f(x)
# based on the analysis of its graphical features.

# The parameters a, b, c, d in the equation correspond to the
# specific points labeled on the x-axis in the provided graph.

# Numerator construction based on the roots of the function.
# Roots are at -b, b, and a double root at d.
numerator = "(x^2 - b^2) * (x - d)^2"

# Denominator construction based on the vertical asymptotes.
# Asymptotes are at x = a (odd multiplicity) and x = c (even multiplicity).
denominator = "(x - a) * (x - c)^2"

# The final equation for f(x) is the ratio of the numerator and the denominator.
# We use '^' to denote exponents for standard mathematical notation.
equation = f"f(x) = ({numerator}) / ({denominator})"

# Print the resulting equation.
print(equation)