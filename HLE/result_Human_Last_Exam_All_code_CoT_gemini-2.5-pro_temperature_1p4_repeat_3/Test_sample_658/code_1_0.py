import math

# The problem asks for the asymptotic value A(B, delta)
# such that d_{B,delta} = Theta(A(B, delta)).
# Based on the analysis, the following expression provides the correct
# asymptotic behavior.

# L is defined as log(delta^-1) in the problem description.
# The formula is B + L / log(L/B + 2).
# Here, we print this formula as a string.

final_expression = "B + L / log(L/B + 2)"
print(final_expression)
