import sympy

# Define the symbols used in the model parameters
phi = sympy.Symbol('phi')
mu = sympy.Symbol('mu')
rho = sympy.Symbol('rho')

# Construct the expression for kappa based on the derivation from the model equations
kappa = - (phi * mu + rho) / (phi * (1 - mu))

# Print the definition of kappa in a readable format
print("The parameter kappa is a composite of other model parameters.")
print("Its definition is:")
sympy.pprint(sympy.Eq(sympy.Symbol('kappa'), kappa))

# Let's break down the final equation for c* as requested
# c* = (kappa * S - 1) / (S - 1)
# The symbols (or 'numbers') in this final equation are c*, kappa, and S.
print("\nIn the final equation for the critical correlation c*:")
print("c* = (kappa * S - 1) / (S - 1)")
print("The 'numbers' or more accurately, the symbols, are:")
print("c*: The critical correlation.")
print("kappa: The parameter defined above.")
print("S: A constant representing the total weighted sum of synaptic proximities.")
