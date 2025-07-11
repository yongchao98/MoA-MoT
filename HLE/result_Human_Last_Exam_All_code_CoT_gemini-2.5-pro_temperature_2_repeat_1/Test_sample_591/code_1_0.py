import sympy

# The purpose of this script is to display the derived definition of kappa (κ)
# using symbolic mathematics for a clear representation.

# Define the symbols for the parameters in the model.
mu, phi, rho, sigma = sympy.symbols('μ φ ρ σ', real=True, positive=True)

# Based on the derivation from the model's steady-state condition in the
# rate-based limit, kappa is defined as the following combination of parameters.
kappa = -(phi * mu**2 + rho * mu) / (phi * sigma**2)

# The final formula for kappa can also be written by factoring out mu.
kappa_factored = -mu * (phi * mu + rho) / (phi * sigma**2)

# Print the final definition of kappa.
print("The definition of κ is:")
sympy.pprint(kappa_factored)

# The final equation includes symbolic variables that were part of the problem description.
# In order to properly display the formula for kappa, it is helpful to use
# greek letters and mathematical symbols. The expression below shows what
# each of these variables are:
# μ (mu): mean firing rate for all synapses
# φ (phi): scaling constant for the presynaptic accumulator
# ρ (rho): offset constant in the Hebbian learning equation
# σ (sigma): standard deviation of synaptic activation (σ² is the variance)