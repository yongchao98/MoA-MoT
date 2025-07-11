import sympy

# Define the symbols used in the model
mu, rho, phi, tau_u, tau_v = sympy.symbols('mu rho phi tau_u tau_v')

# Define the expression for kappa based on the derivation
kappa_expression = - (phi * mu + rho) / phi * (tau_u + tau_v) / (tau_u * tau_v)

# We can also write it in a different form
kappa_alternative_form = - (mu + rho / phi) * (1 / tau_u + 1 / tau_v)

# Print the definition of kappa
print("The definition of the parameter kappa (κ) is:")
print(f"κ = {sympy.pretty(kappa_alternative_form, use_unicode=True)}")

# Explain the parameters involved in the definition
print("\nWhere the terms in the equation are:")
print("μ (mu)   : The mean firing rate of presynaptic neurons.")
print("ρ (rho)  : The offset constant in the Hebbian learning rule.")
print("φ (phi)  : The scaling constant for the presynaptic accumulator.")
print("τ_u (tau_u): The time constant for the postsynaptic accumulator.")
print("τ_v (tau_v): The time constant for the presynaptic accumulator.")
