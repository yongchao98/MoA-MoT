import sys

# The Born rule is a fundamental postulate of quantum mechanics that gives the probability
# that a measurement of a quantum system will yield a given result. In its simplest form,
# it states that the probability density of finding a particle at a given point, when measured,
# is proportional to the square of the magnitude of the particle's wavefunction at that point.
# A key consequence is the normalization condition, which states the total probability is 1.

# Using Unicode characters for mathematical symbols:
# ∫ (Integral): \u222b
# ψ (psi): \u03c8
# ² (Superscript Two): \u00b2
# τ (tau, for a volume element): \u03c4

name = "Max Born"
# The equation represents the normalization condition for a wave function.
# It states that the integral of the squared magnitude of the wave function over all space is equal to 1.
# This ensures that the particle exists somewhere.
integral = "\u222b"
psi_squared = "|\u03c8|\u00b2"
d_tau = "d\u03c4"
equals = "="
one = "1"

# We print the name and then construct and print the equation part by part.
# The user wants each number in the final equation printed.
print(f"Name: {name}")
sys.stdout.write("Equation: ")
sys.stdout.write(f"{integral} ")
sys.stdout.write(f"{psi_squared} ")
sys.stdout.write(f"{d_tau} ")
sys.stdout.write(f"{equals} ")
sys.stdout.write(f"{one}\n")