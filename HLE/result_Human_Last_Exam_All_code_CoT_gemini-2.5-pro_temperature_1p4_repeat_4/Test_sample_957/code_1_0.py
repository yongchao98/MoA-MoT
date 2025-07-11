import sympy as sp

# Define the symbols
rho_m, tau_m, k0, d = sp.symbols('rho_m tau_m k_0 d', complex=True)
I = sp.I

# Phase term for a single pass
phi_d = k0 * d
# Phase term for a round trip
phi_2d = 2 * k0 * d

# Denominator term (from geometric series summation)
denominator = 1 - rho_m**2 * sp.exp(I * phi_2d)

# Overall Transmission Coefficient (tau)
tau_numerator = tau_m**2 * sp.exp(I * phi_d)
tau = tau_numerator / denominator

# Overall Reflection Coefficient (rho)
rho_numerator = rho_m - rho_m * (rho_m**2 - tau_m**2) * sp.exp(I * phi_2d)
rho = rho_numerator / denominator

# Print the final expressions
print("Derived expressions for the Fabry-Pérot cavity coefficients:\n")
print("Transmission coefficient (tau):")
sp.pretty_print(tau)
print("\nReflection coefficient (rho):")
sp.pretty_print(rho)

print("\nComparing with the answer choices:")
print("The derived expression for tau matches option D.")
print("The derived expression for rho does not exactly match any option, likely due to typos in the choices.")
print("The expression for tau in option D is standard and correct, making D the most likely intended answer.\n")

# Print the formulas from option D for comparison
print("Option D provides:")
tau_D_str = "τ = (τ_m**2 * exp(i*k_0*d)) / (1 - ρ_m**2 * exp(i*2*k_0*d))"
rho_D_str = "ρ = (1 - (ρ_m - τ_m**2) * exp(i*2*k_0*d) * ρ_m) / (1 - ρ_m**2 * exp(i*2*k_0*d))"
print(tau_D_str)
print(rho_D_str)
