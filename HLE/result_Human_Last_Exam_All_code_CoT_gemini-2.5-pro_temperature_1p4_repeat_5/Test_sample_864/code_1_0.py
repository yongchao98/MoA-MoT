import sympy

# Define the D-transform variable
D = sympy.Symbol('D')

# Based on the analysis, a likely intended simplified model for the
# causal spectral factor is G(D) = 1 + 2*D.
# The corresponding anti-causal factor is G(D^-1).
# Let's define G_inv_D which corresponds to G(D^{-1})
G_inv_D = 1 + 2 * D**-1

# The whitening filter W(D) that makes the resulting channel causal is 1/G(D^{-1}).
W_D = 1 / G_inv_D

# Simplify the expression for the whitening filter W(D)
W_D_simplified = sympy.simplify(W_D)

# The resulting filter is W(D)
# We can print the components of the equation for clarity
# The final answer is the expression for W(D)
# The problem asks for the filter W(D)
# W(D) = 1 / G(D^{-1}) where G(D^{-1}) = 1 + 2*D^{-1}
# So W(D) = 1 / (1 + 2*D^{-1})

# To show the output clearly, we'll represent the filter components
numerator, denominator = sympy.fraction(W_D_simplified)

print("The whitening filter W(D) is found by taking the inverse of the anti-causal spectral factor G(D^-1).")
print("Assuming a simplified channel model where G(D) = 1 + 2*D, then G(D^-1) = 1 + 2*D^-1.")
print("So, the whitening filter is W(D) = 1 / (1 + 2*D^-1)")
print("This can be simplified to:")
print(f"W(D) = {sympy.pretty(W_D_simplified)}")

final_expression = f"1/({sympy.pretty(G_inv_D)})"

# Let's break down the final equation as requested.
# W(D) = 1 / G(D^{-1})
# Our G(D^{-1}) is 1 + 2*D^{-1}
# So we print the numerator and denominator of this expression.
num = 1
g0 = 1
g1 = 2

print("\nThe final equation is of the form W(D) = Numerator / Denominator")
print(f"Numerator = {num}")
print(f"Denominator = g0 + g1*D^-1")
print(f"g0 = {g0}")
print(f"g1 = {g1}")
print(f"Thus, the equation is: {num} / ({g0} + {g1}*D^-1)")
