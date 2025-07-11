import sympy

# Define the symbolic variables used in the equation
x = sympy.Symbol('x')  # End-to-end separation
l = sympy.Symbol('l', positive=True)  # Length of a single link
n = sympy.Symbol('n', integer=True, positive=True) # Number of links
E0 = sympy.Symbol('E(0)') # Kinetic energy at zero extension

# Based on the derivation, the relationship between energy E and extension x is:
# E(x) = E(0) * exp(x**2 / (n**2 * l**2))
# The force is F = -dE/dx.

# Construct the final equation for the force F(x):
# F(x) = - (2 * E(0) * x) / (n**2 * l**2) * exp(x**2 / (n**2 * l**2))

coefficient = -2
main_term = (E0 * x) / (n**2 * l**2)
exponent_term = (x**2) / (n**2 * l**2)

force_expression = coefficient * main_term * sympy.exp(exponent_term)

print("The derived force law F(x) for a thermally isolated freely jointed chain is:")
print("-" * 70)

# Print the equation in a clear, formatted string.
# The prompt is specific: "output each number in the final equation!"
print(f"F(x) = -(2 * E(0) * x / (n**2 * l**2)) * exp(x**2 / (n**2 * l**2))")

print("\nLet's break down the components of this equation to be explicit about the numbers:")
print("The overall numerical coefficient is -2.")
print("The term 'n**2 * l**2' appears in both the pre-factor and the exponent. The number in this term is the exponent 2.")
print("The term 'x**2' in the exponent also involves the number 2 as the exponent.")

print("\nFor a cleaner mathematical representation, here is the output from the sympy library:")
final_eq = sympy.Eq(sympy.Symbol("F(x)"), force_expression)
sympy.pprint(final_eq, use_unicode=True)
print("-" * 70)
