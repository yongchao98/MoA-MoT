import sys
import io

# Set stdout to utf-8 to ensure special characters are printed correctly.
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# The user wants to complete a list of physicists and their equations.
# The entry to be completed is for "Max Born".
# Max Born's most significant contribution to quantum mechanics is the Born rule,
# which provides the probabilistic interpretation of the wave function.
# The rule states that the probability density (P) of finding a particle at a
# position (x) is equal to the square of the absolute value of its wave function (ψ).
# The equation is P(x) = |ψ(x)|².

physicist_name = "Max Born"

# The number in the equation is 2 (the exponent).
exponent = 2

# We format the equation in LaTeX style, similar to the user's examples.
equation_string = f"\(P(x) = |\\psi(x)|^{exponent}\)"

# Print the name and the equation.
print(f"Name: {physicist_name}")
print(f"Equation: {equation_string}")