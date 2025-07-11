import sympy

# Define the symbols
mu, mu_0, N, w, g, x1, x2, I1, I2 = sympy.symbols('μ μ_0 N w g x_1 x_2 I_1 I_2')

# Construct the expression for the work done by the source
# Based on the derivation: W = -1/2 * (L(x2) - L(x1)) * (I2**2 - I1**2)
# where L(x2) - L(x1) = (N**2 * w / g) * (mu - mu_0) * (x2 - x1)

L_diff = (N**2 * w / g) * (mu - mu_0) * (x2 - x1)
I_sq_diff = I2**2 - I1**2
W = -sympy.Rational(1, 2) * L_diff * I_sq_diff

# Format the final equation to be readable
final_equation_str = f"W = - (({mu} - {mu_0}) / (2*{g})) * {N}**2 * {w} * ({x2} - {x1}) * ({I2}**2 - {I1}**2)"

print("The final equation for the work W is:")
print(final_equation_str)

# To satisfy the requirement of outputting each 'number' (component) in the equation:
print("\nThe components of the final equation are:")
print(f"Overall sign: -1")
print(f"Permeability term: ({mu} - {mu_0})")
print(f"Geometric term in denominator: 2*{g}")
print(f"Coil and geometry term in numerator: {N}**2 * {w}")
print(f"Displacement term: ({x2} - {x1})")
print(f"Current term: ({I2}**2 - {I1}**2)")