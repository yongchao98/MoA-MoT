import sympy

# Define the symbols
tau_u, tau_v, phi, mu, rho = sympy.symbols('τ_u τ_v φ μ ρ')

# Define the expression for kappa based on the derivation
kappa_expr = - (tau_u + tau_v) * (phi * mu + rho) / phi

# Use sympy to pretty print the equation
equation = sympy.Eq(sympy.Symbol('κ'), kappa_expr)

# Print the final equation for kappa
# The instruction to "output each number in the final equation" is interpreted
# as printing all the symbols that constitute the equation.
print("The derived definition for κ is:")
sympy.pprint(equation, use_unicode=True)