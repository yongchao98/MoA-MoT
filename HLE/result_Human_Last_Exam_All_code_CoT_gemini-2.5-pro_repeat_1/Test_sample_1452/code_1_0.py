import sys

# This script will print the famous equation associated with Max Born.
# The Born rule is a fundamental postulate of quantum mechanics which gives the probability
# that a measurement of a quantum system will yield a given result.
# It states that the probability density of finding a particle at a given point is
# proportional to the square of the magnitude of the particle's wavefunction at that point.

# The wavefunction is represented by the symbol 'ψ' (psi).
# The equation is the probability density P = |ψ|²
# The '2' in the equation is the exponent.

name = "Max Born"
# Using a raw string to handle the backslash for psi
equation_text = "|\\psi|^2"

print(f"Name: {name}")
# We need to output the final equation, including the number '2'.
print(f"Equation: {equation_text}")