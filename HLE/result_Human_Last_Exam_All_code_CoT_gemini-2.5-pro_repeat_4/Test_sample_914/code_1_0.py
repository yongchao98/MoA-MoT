import sympy

# Define the symbolic variables for the system parameters
a, D, s, mu_0, I_0, sigma_1, sigma_2 = sympy.symbols('a D s mu_0 I_0 sigma_1 sigma_2', real=True, positive=True)

# The x-directed total force on the conducting material as given in Answer A
Fx = -a * D * (mu_0 / 2) * (I_0**2 / D**2) * (sigma_2 / (sigma_1 + sigma_2))**2

# Print the formula for the force
print("The x-directed total force on the conducting material is given by the equation:")
sympy.pprint(Fx, use_unicode=True)

# To display it in a more standard way:
print("\nOr, written more conventionally:")
print(f"F_x = -a*D * (μ₀/2) * (I₀² / D²) * (σ₂ / (σ₁ + σ₂))²")
