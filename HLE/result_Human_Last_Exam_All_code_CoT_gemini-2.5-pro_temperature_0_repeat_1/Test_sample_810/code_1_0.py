import sympy

# Define the symbols
c, t = sympy.symbols('c t')
K = sympy.Function('K')(sympy.Function('gamma')(t))
theta = sympy.Function('theta')(t)

# The derived expression for the derivative of theta(t)
theta_prime_expr = c * sympy.cos(theta)**2 + (1/c) * K * sympy.sin(theta)**2

# Print the final equation
# The problem asks to output each number in the final equation.
# The numbers are the powers '2' and the coefficients 'c' and '1/c'.
print("The derivative of theta(t) is given by the equation:")
sympy.pprint(sympy.Eq(sympy.Derivative(theta, t), theta_prime_expr), use_unicode=False)

print("\nBreaking down the equation:")
print(f"The first term is: {c} * cos^2(theta(t))")
print(f"The second term is: (1/{c}) * K(gamma(t)) * sin^2(theta(t))")