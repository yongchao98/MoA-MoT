import sympy

# Define the symbols based on the problem description
phi = sympy.Symbol('phi')
mu = sympy.Symbol('mu')
rho = sympy.Symbol('rho')
tau_u = sympy.Symbol('tau_u')
tau_v = sympy.Symbol('tau_v')

# The derived expression for kappa
kappa_expr = - (phi * mu + rho) * (tau_u + tau_v) / phi

# Print the definition of kappa as an equation
print(sympy.Eq(sympy.Symbol('kappa'), kappa_expr))