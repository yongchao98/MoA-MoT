import sympy as sp

# Define the symbols
mu_0, K_0, omega, t, omega_p, d, c = sp.symbols('mu_0 K_0 omega t omega_p d c', real=True, positive=True)
i_x = sp.Symbol('i_x') # Represents the unit vector in x

# Construct the expression for the force per unit area based on option E
# The core derived part is mu_0 * K_0**2 * cos(omega*t)**2 / (2 * cosh(omega_p*d/c)**2)
# Option E includes an additional exponential term
force_expression = (sp.Rational(1, 2)) * (mu_0 * K_0**2 * sp.cos(omega*t)**2 / (sp.cosh(omega_p * d / c)**2)) * sp.exp(-omega * d / c)

# Print the final formatted equation
print("The force per unit area on the x = d plane is:")
print(f"f_vec = i_x * ({force_expression})")
print("\nBreaking down the expression from option E:")
print(f"Direction: i_x")
print(f"Coefficient: 1/2")
print(f"Magnetic permeability term: mu_0")
print(f"Surface current term: K_0^2")
print(f"Time-dependent term: cos^2(omega*t)")
print(f"Denominator term: cosh^2(omega_p*d/c)")
print(f"Exponential perturbation term: exp(-omega*d/c)")
