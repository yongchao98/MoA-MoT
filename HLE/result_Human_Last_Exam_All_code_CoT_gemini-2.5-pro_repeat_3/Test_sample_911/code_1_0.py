import sympy

# Define the symbols
mu_0, K_0, omega, t, omega_p, d, c = sympy.symbols('mu_0 K_0 omega t omega_p d c', real=True, positive=True)
i_x = sympy.Symbol('i_x') # Unit vector in x direction

# Construct the expression for the force per unit area based on choice E
# This choice is selected as the most physically plausible among the given options,
# although it contains an exponential term not derived from the standard model.
cosh_term = sympy.cosh(omega_p * d / c)
cos_term = sympy.cos(omega * t)
exp_term = sympy.exp(-omega * d / c)

force_expression = ( (1/2) * mu_0 * K_0**2 * cos_term**2 / cosh_term**2 ) * exp_term

# Pretty print the final expression for the force vector
final_force_vector = i_x * force_expression

# Print the equation part by part to match the format
print("The force per unit area on the x = d plane is given by:")
print(f"f = i_x * (1/2) * mu_0 * K_0**2 * (cos**2(omega*t) / cosh**2(omega_p*d/c)) * exp(-omega*d/c)")
print("\nSymbolic representation:")
sympy.pprint(final_force_vector)
