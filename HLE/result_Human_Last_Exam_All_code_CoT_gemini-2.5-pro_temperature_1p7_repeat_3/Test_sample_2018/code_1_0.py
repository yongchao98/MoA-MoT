import math

# The problem asks for the value of the computational factor 'C' in the 
# Carman-Kozeny source term used in early simulation-only papers on melting.
# S = -C * ((1-f)^2 / (f^3 + ε)) * u
# where f is the liquid fraction and ε is a small constant.

# Based on the foundational 1988 paper by Brent, Voller, and Reid, this
# factor 'C' was set to a specific large value for their simulation-only work.
C = 1.6e6

# The small constant in the denominator for numerical stability.
# The problem states this is on the order of 10^-3.
epsilon = 1e-3

# Let's represent the variables in the equation with symbolic placeholders for clarity.
f = 'f'  # Liquid fraction (porosity)
u = 'u'  # Velocity vector

print("The Carman-Kozeny source term is often expressed as:")
print(f"S = -C * (1 - {f})^2 / ({f}^3 + ε) * {u}\n")
print("According to the prior published implementation, the values used were:")

print(f"Computational Factor (C): {C:.1e}")
print(f"Denominator Constant (ε): {epsilon}")

print("\nSubstituting the constant values into the equation gives:")
print(f"S = -{C:.1e} * (1 - {f})^2 / ({f}^3 + {epsilon}) * {u}")
