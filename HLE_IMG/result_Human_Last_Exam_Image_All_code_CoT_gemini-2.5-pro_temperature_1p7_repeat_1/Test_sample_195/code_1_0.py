# Based on the graphical analysis, we determine the components of the rational function.
# The roots of the function are at x = d, x = b, and x = -b.
# The vertical asymptotes are at x = a and x = c.
# These correspond to the factors of the numerator and denominator, respectively.
# The problem asks for the equation, including all the labeled points.

# We represent the parameters with their character labels as in the graph.
param_a = 'a'
param_b = 'b'
param_c = 'c'
param_d = 'd'

# We construct the equation string using f-string formatting.
# The numerator is (x-d)(x-b)(x+b)
# The denominator is (x-a)(x-c)
# The expression is formatted for clarity.
equation = f"f(x) = ((x - {param_d})*(x - {param_b})*(x + {param_b})) / ((x - {param_a})*(x - {param_c}))"

print(equation)