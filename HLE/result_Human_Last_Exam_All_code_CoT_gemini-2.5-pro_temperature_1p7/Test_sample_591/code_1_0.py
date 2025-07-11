# The user asks for the definition of kappa in the given context.
# Based on the derivation from the model's equations, kappa is a composite
# parameter that depends on several other parameters of the model.

# This script prints the derived definition of kappa as a string.
kappa_definition = "kappa = -(tau_u + tau_v) * (mu + rho / phi)"

print("Based on a steady-state analysis of the provided model equations, the definition of kappa is:")
print(kappa_definition)