import sympy

# This Python code symbolically derives and prints the equation for the shape of the set S
# under the simplifying assumption that the vectors {y_i} are orthogonal.

# Let's consider n=3 for a clear, generalizable example.
n = 3

# Symbolic representation for the components of a vector in S
x_symbols = sympy.symbols(f'x_1:{n+1}')

print("Based on the analysis, if the vectors {y_i} are orthogonal, the set S forms a simplex.")
print("A simplex is a surface defined by a single linear equation, restricted to the non-negative orthant.")
print("\nThe equation for this simplex is derived as follows:")

# Build the equation string
equation_parts = []
for i in range(n):
    # Add the term x_i / ||y_i||^2
    equation_parts.append(f"x_{i+1} / ||y_{i+1}||^2")
    # Add a plus sign if it's not the last term
    if i < n - 1:
        equation_parts.append("+")

# Join the parts and add the equality
final_equation = " ".join(equation_parts) + " = 1"

print(final_equation)

print("\nTo satisfy the output format of printing each number/symbol, here is the equation again, broken down:")
# "each number" is interpreted as each symbolic part of the equation.
# Since the equation is symbolic, we print the symbols themselves.
equation_components = []
for i in range(n):
    equation_components.append(f"x_{i+1}")
    equation_components.append("/")
    equation_components.append(f"||y_{i+1}||^2")
    if i < n - 1:
        equation_components.append("+")
equation_components.append("=")
equation_components.append("1")

for part in equation_components:
    print(part)