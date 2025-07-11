import math

# This script presents the final derived equation for the amplitude of the electric field.

# Define symbolic representations of the physical quantities
A = "A"  # Initial amplitude of the incident wave
alpha = "α" # The rate of change of the slab's electromagnetic properties
L = "L"  # The length of the slab in the x-direction
c = "c"  # The speed of light in vacuum
e = "e"  # Euler's number, the base of the natural logarithm

# The derived amplitude of the electric field at the rightmost boundary (x=L)
# is denoted by A_L. The equation shows its dependence on the initial
# amplitude and the properties of the slab.

print("The final equation for the amplitude of the electric field (A_L) at the boundary x=L is:")

# The final equation is A_L = A * exp(-(alpha * L) / c)
# We print each part of the equation clearly.
final_equation = f"A_L = {A} * {e}^(-({alpha} * {L}) / {c})"
print(final_equation)

print("\nWhere the terms in the equation are:")
print(f" A_L: The final amplitude at x=L.")
print(f"  {A}: The initial amplitude of the wave before entering the slab.")
print(f"  {e}: The base of the natural logarithm (approx. {math.e}).")
print(f"  {alpha}: The time-variation coefficient from the problem statement.")
print(f"  {L}: The length of the slab.")
print(f"  {c}: The speed of light in vacuum.")
print("  -1: The numeric constant in the exponent, indicating amplitude decay for α > 0.")