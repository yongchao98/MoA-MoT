import sys
import io

# Set stdout to utf-8 to handle special characters
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# The final equation for Max Born is the Born rule, which gives the probability density.
# The probability density (ρ) is equal to the squared magnitude of the wave function (Ψ).
# This is a fundamental principle of quantum mechanics.

name = "Max Born"
# Using ρ to represent probability density and Ψ for the wave function.
# The number in the equation is 2 (the exponent).
equation_latex = "\\rho = |\\Psi|^2"

print(f'Name: "{name}"')
# The user wants the equation to be displayed with LaTeX formatting
print(f'Equation: \\({equation_latex}\\)')
