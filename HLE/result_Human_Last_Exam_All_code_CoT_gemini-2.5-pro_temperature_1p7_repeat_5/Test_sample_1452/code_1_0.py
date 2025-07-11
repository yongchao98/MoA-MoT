import sys
# The user wants to find the famous equation associated with the physicist Max Born.
# Max Born's most fundamental contribution to quantum mechanics is the Born rule.
# The rule provides a way to calculate the probability of finding a quantum particle
# in a certain state. It states that the probability density is equal to the
# square of the magnitude of the wave function, denoted as |Ψ|².
# This equation includes the number '2' as the exponent.
# The following code prints the name and the corresponding equation in LaTeX format.

name = "Max Born"
# The equation for the Born rule. We use a raw string (r"...") for the LaTeX format
# to avoid issues with backslashes.
equation = r"\( |\Psi|^2 \)"

print(f"Name: {name}")
print(f"Equation: {equation}")