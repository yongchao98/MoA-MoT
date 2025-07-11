import sympy

# Define the symbols
x, l, n = sympy.symbols('x ell n')
E0 = sympy.Symbol('E(0)')

# Construct the force law expression based on the derivation
# The force is attractive (negative sign).
# It's proportional to the initial kinetic energy E(0) and the extension x.
# The stiffness depends on the number of links n and the segment length l.
# There is a non-linear stiffening term from the exponential.
coefficient = - (3 * E0) / (n**2 * l**2)
linear_term = x
exponential_term = sympy.exp((3 * x**2) / (2 * n**2 * l**2))

force_law = coefficient * linear_term * exponential_term

# Print the formula for the force law
print("The force law F(x) is:")
sympy.pprint(force_law, use_unicode=True)

# For the final answer format, we print the equation explicitly.
# The instruction "output each number in the final equation" is interpreted as
# printing the components of the derived formula.
print("\nFinal Equation:")
print(f"F(x) = - (3 * E(0) / (n**2 * l**2)) * x * exp(3 * x**2 / (2 * n**2 * l**2))")
final_answer_string = f"- (3 * E(0) / (n**2 * l**2)) * x * exp(3 * x**2 / (2 * n**2 * l**2))"
print(f"<<<{final_answer_string}>>>")
