import sympy

# Define the symbols
tau, rho, tau_m, rho_m, k0, d = sympy.symbols('τ ρ τ_m ρ_m k_0 d')
i = sympy.I

# Expression for the overall transmission coefficient
tau_expr = (tau_m**2 * sympy.exp(i * k0 * d)) / (1 - rho_m**2 * sympy.exp(i * 2 * k0 * d))

# Expression for the overall reflection coefficient
rho_expr = (1 - (rho_m - tau_m**2) * sympy.exp(i * 2 * k0 * d) * rho_m) / (1 - rho_m**2 * sympy.exp(i * 2 * k0 * d))

# Print the equations
print("The overall transmission coefficient is:")
print(f"τ = {sympy.pretty(tau_expr, use_unicode=False)}")
print("\nThe overall reflection coefficient is:")
# The pretty print for the rho expression is a bit wide, let's just print the normal string representation
print(f"ρ = (1 - (ρ_m - τ_m**2)*exp(I*2*k_0*d)*ρ_m) / (1 - ρ_m**2*exp(I*2*k_0*d))")