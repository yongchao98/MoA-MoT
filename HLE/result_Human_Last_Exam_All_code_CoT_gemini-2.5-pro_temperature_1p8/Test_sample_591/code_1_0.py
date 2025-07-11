import sympy

# Define the symbolic variables from the model
tau_u, tau_v, mu, rho, phi = sympy.symbols('tau_u tau_v mu rho phi')

# The derived expression for kappa
kappa_expression = -(tau_u + tau_v) * (mu + rho / phi)

# Print the result in a readable format
print("The definition of kappa is:")
sympy.pprint(kappa_expression, use_unicode=False)

# To express the final answer in the requested format, we format the string.
final_answer_str = str(kappa_expression).replace(' ', '')
print(f"\nFinal Answer String: {final_answer_str}")