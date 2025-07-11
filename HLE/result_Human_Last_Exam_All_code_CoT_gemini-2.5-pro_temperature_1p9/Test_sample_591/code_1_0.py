import sympy

# This script will define the variable kappa based on the provided biophysical model.
# The definition is derived algebraically from the steady-state conditions of the system.
# We will print the final equation and explain the meaning of each parameter.

# Define the symbols using sympy for clear mathematical representation.
# Using unicode characters for better readability.
tau_u, tau_v, mu, rho, phi, kappa = sympy.symbols('τ_u, τ_v, μ, ρ, φ, κ')

# The derived expression for kappa is:
kappa_definition = - (tau_u + tau_v) * (mu + rho / phi)

# Create a formatted string for the equation.
# Using sympy.pretty() for a nicely formatted output.
equation_str = sympy.pretty(sympy.Eq(kappa, kappa_definition), use_unicode=True)

# Print the final result
print("Based on the steady-state analysis of the dendritic plasticity model, the definition of the parameter κ is:")
print()
print(equation_str)
print()
print("Where the parameters are:")
print(f"  {tau_u}: Time constant of the postsynaptic accumulator")
print(f"  {tau_v}: Time constant of the presynaptic accumulator")
print(f"  {mu}: Mean firing rate of presynaptic neurons")
print(f"  {rho}: Offset constant in the synaptic efficacy equation")
print(f"  {phi}: Scaling constant in the presynaptic accumulator equation")