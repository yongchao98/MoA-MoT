import sympy as sp

# Define symbols
A, B, theta = sp.symbols('A B theta')

# From the derivation, the corrector term is exp(psi_1(theta))
# psi_1'(theta) = A*tan(theta/2) + B
# Integrating gives psi_1(theta) = -2*A*ln(cos(theta/2)) + B*theta
psi_1 = -2 * A * sp.log(sp.cos(theta / 2)) + B * theta

# The corrector is exp(psi_1)
corrector = sp.exp(psi_1)
corrector_simplified = sp.exp(B * theta) * (sp.cos(theta / 2))**(-2 * A)

# Print the result in a readable format
print("The corrector to the large-distance behavior of omega is:")
print(sp.pretty(corrector_simplified, use_unicode=True))

# The final asymptotic behavior of omega(r,theta) is proportional to:
# 1/sqrt(r) * exp(-r*(1-cos(theta))) * corrector
# We can print the components of the corrector's power expression.
# The expression is of the form (cos(theta/2))^(p1) * exp(p2)
p1 = -2 * A
p2 = B * theta

print("\nThe structure is (cos(theta/2)) to the power of P1, multiplied by exp(P2).")
print(f"P1 = {p1}")
print(f"P2 = {p2}")
