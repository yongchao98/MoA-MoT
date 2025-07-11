import sympy as sp

# Define the symbols
mu_0, K_0, omega, t, omega_p, d, c = sp.symbols('mu_0 K_0 omega t omega_p d c', real=True, positive=True)
i_x = sp.Matrix([1, 0, 0])

# Construct the expression for the force per unit area from option E
cosh_arg = omega_p * d / c
exp_arg = -omega * d / c
force_magnitude = (sp.Rational(1, 2) * mu_0 * K_0**2 * sp.cos(omega*t)**2) / (sp.cosh(cosh_arg)**2) * sp.exp(exp_arg)

# Express the force as a vector in the x-direction
f_vec = force_magnitude * i_x.T

# Print the formula in a readable format
# For display, we use sympy's pretty print, but the output for the user will be a string
force_expression_str = f"f = i_x * (1/2) * (mu_0 * K_0**2 * cos(omega*t)**2) / (cosh({sp.pretty(cosh_arg)})**2) * exp({sp.pretty(exp_arg)})"

print("The force per unit area on the x=d plane is given by the expression from option E:")
# This part of the code aims to display the final answer in a structured way as requested by the user prompt
print(f"   f_vec = i_x * (1/2) * (mu_0 * K_0**2 * cos(omega*t)**2) / (cosh(omega_p*d/c)**2) * exp(-omega*d/c)")
print("\nLet's print each component of the final equation to be clear:\n")
print(f"Direction vector: i_x")
print(f"Factor: 1/2")
print(f"Magnetic permeability constant: mu_0")
print(f"Surface current amplitude squared: K_0**2")
print(f"Time-dependent term: cos^2(omega*t)")
print(f"Denominator term: cosh^2(omega_p*d/c)")
print(f"Exponential term: exp(-omega*d/c)")