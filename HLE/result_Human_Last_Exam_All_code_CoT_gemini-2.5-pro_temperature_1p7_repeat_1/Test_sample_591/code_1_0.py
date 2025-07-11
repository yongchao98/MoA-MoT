import sympy

# Define the symbols used in the model
rho, phi, mu, tau_u, tau_v = sympy.symbols('rho phi mu tau_u tau_v')

# The derived definition for kappa
# kappa represents a combination of biophysical parameters that determines
# the requirement for correlated activity to stabilize synapses.
kappa_expression = -(rho / phi + mu) * (tau_u + tau_v)

# Print the definition of kappa
# The instruction "output each number in the final equation" is interpreted as
# outputting all the variables involved in the definition.
print("The definition of kappa is:")
print(f"kappa = {kappa_expression}")
