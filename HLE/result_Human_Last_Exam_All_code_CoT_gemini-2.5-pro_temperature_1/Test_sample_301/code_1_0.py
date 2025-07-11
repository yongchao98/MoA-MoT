import sympy

# Define the symbols for the concentrations and rate constants
B, F = sympy.symbols('[B] [F]')
k2, k3, k4, k5 = sympy.symbols('k2 k3 k4 k5')

# The relationship is a direct proportionality along the pathway
# B -> Malonyl-CoA -> Acetyl-CoA -> Pyruvate -> F
proportionality_expression = k2 * k3 * k4 * k5

# Print the final derived relationship in a clear format
print("The relationship between the concentration of PEP ([F]) and 3-Hydroxypropionate ([B]) is determined by the direct pathway connecting them.")
print("The pathway is: 3-Hydroxypropionate -k2-> Malonyl-CoA -k3-> Acetyl-CoA -k4-> Pyruvate -k5-> PEP")
print("\nAssuming a direct proportionality for each step, the final expression is derived by multiplying the rate constants along this path.")
print("\nThe final relationship is:")
# Using sympy.pretty_print for a nicer output
# We create an equation to show the proportionality
final_equation = sympy.Eq(F, B * proportionality_expression)

# We want to print F is proportional to B * k2 * k3 * k4 * k5
# Let's format it manually for clarity
print(f"{F} ∝ {B} * {k2} * {k3} * {k4} * {k5}")