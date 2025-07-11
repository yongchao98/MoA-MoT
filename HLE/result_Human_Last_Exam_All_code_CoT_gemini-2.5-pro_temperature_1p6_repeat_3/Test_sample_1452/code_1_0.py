import sys

# Name of the physicist
name = "Max Born"

# Max Born's most famous contribution is the Born rule, which provides the
# probabilistic interpretation of the wave function in quantum mechanics.
# The equation states that the probability density (P) of finding a particle
# at a position (x) is equal to the square of the absolute value of the
# wave function psi(x).
# The number '2' is an essential part of this equation.
equation_latex = r"P(x) = |\psi(x)|^2"

# Print the name and the equation in the specified format
print(f"Name: {name}")
print(f"Equation: \\( {equation_latex} \\)")
