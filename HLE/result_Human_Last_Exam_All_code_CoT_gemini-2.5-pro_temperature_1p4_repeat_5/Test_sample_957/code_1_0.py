import sympy

# Define the symbols used in the Fabry-Pérot equations
tau_m, rho_m, k0, d = sympy.symbols('tau_m rho_m k_0 d', real=True)
i = sympy.I

# Expression for the overall transmission coefficient (tau)
# This is derived by summing the geometric series of transmitted partial waves.
denominator = 1 - rho_m**2 * sympy.exp(i * 2 * k0 * d)
tau_numerator = tau_m**2 * sympy.exp(i * k0 * d)
tau = tau_numerator / denominator

# Expression for the overall reflection coefficient (rho)
# This is derived by summing the directly reflected wave and the series of waves that exit back through the first mirror.
rho_numerator = rho_m - (rho_m**2 - tau_m**2) * rho_m * sympy.exp(i * 2 * k0 * d)
rho = rho_numerator / denominator

# Print the final derived equations
print("Derived Transmission Coefficient (τ):")
# The following line uses sympy's pretty print for a clear mathematical representation
# For plain text, one could use str(tau)
sympy.pprint(tau, use_unicode=True)
print("\n")
print("Derived Reflection Coefficient (ρ):")
sympy.pprint(rho, use_unicode=True)

# As text for easier reading of the final result:
print("\n--- Simplified Formulas ---")
print(f"τ = (τ_m**2 * exp(i*k_0*d)) / (1 - ρ_m**2 * exp(i*2*k_0*d))")
print(f"ρ = (ρ_m - (ρ_m**2 - τ_m**2)*ρ_m*exp(i*2*k_0*d)) / (1 - ρ_m**2 * exp(i*2*k_0*d))")
print("\nNote: These derived formulas correspond to option D after correcting likely typos in its expression for ρ.")
