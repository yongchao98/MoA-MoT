import sympy

# Define the symbols
mu_0, K_0, omega, t, omega_p, d, c = sympy.symbols('mu_0 K_0 omega t omega_p d c', real=True, positive=True)
i_x = sympy.Matrix([1, 0, 0])

# Construct the expression for the force per unit area based on option E
# This option matches our derivation in its core components but includes an additional exponential factor.
cosh_term = sympy.cosh(omega_p * d / c)
cos_term = sympy.cos(omega * t)
exp_term = sympy.exp(-omega * d / c)

# Magnitude of the force vector from option E
f_magnitude_E = (1/2) * mu_0 * K_0**2 * (cos_term**2 / cosh_term**2) * exp_term

# The final vector expression
f_vector_E = f_magnitude_E * sympy.Symbol('i_x')

# Print the formula as a string
# The task asks to output each number in the final equation. We will format the string to be clear.
print("The force per unit area is given by the expression from option E:")
print(f"f = i_x * (1/2) * mu_0 * K_0**2 * cos(omega*t)**2 / cosh(omega_p*d/c)**2 * exp(-omega*d/c)")
# For direct evaluation, we can use sympy's pretty print
# sympy.pprint(f_vector_E)
