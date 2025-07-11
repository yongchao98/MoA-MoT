import sys

# The physicist's name
name = "Max Born"

# Max Born's famous equation is the Born rule, which gives the probability density
# of finding a particle at a given point. P is the probability density, and Ψ
# is the wave function. The equation contains the number 2.
equation_text = r"P = |\Psi|^2"

# Print the final output in the required format.
# Using LaTeX-style delimiters as in the user's examples.
print(f"Name: {name}")
print(f"Equation: \\( {equation_text} \\)")
