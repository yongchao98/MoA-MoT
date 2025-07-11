import sympy as sp

# This script symbolically represents the derived relationship for pi_0.

# Define the symbols used in the equation
# pi_0 represents the steady-state probability of being in state 0
# rho is the ratio of the arrival rate lambda to the service rate mu
pi_0 = sp.Symbol('π₀')
rho = sp.Symbol('ρ')

# The derived final equation is pi_0 = e**(-rho).
# We represent this using sympy's Eq for equality and exp for the exponential function.
final_equation = sp.Eq(pi_0, sp.exp(-rho))

# Print the final equation in a readable format.
# sp.pprint allows for nice printing of mathematical expressions, including Unicode characters.
# The final equation includes the symbols π₀ and ρ, which represent the "numbers" in this symbolic context.
print("The final equation for the steady-state probability π₀ is:")
sp.pprint(final_equation, use_unicode=True)