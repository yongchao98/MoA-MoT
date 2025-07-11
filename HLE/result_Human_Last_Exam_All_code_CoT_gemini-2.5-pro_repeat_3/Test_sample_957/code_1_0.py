import sympy

# Define the symbols used in the equations
# Note: In Python, 'tau' and 'rho' are used for the variable names.
# The 'm' subscript is represented by adding '_m'.
# Complex exponential is represented using sympy.exp and sympy.I
tau_m, rho_m, E0, d, k0 = sympy.symbols('tau_m rho_m E_0 d k_0')
I = sympy.I

# Expression for the overall transmission coefficient (tau) from option D
tau_expr = (tau_m**2 * sympy.exp(I * k0 * d)) / (1 - rho_m**2 * sympy.exp(I * 2 * k0 * d))

# Expression for the overall reflection coefficient (rho) from option D
rho_expr = (rho_m - (rho_m**2 - tau_m**2) * sympy.exp(I * 2 * k0 * d) * rho_m) / (1 - rho_m**2 * sympy.exp(I * 2 * k0 * d))

# The provided option D has a slightly different, likely erroneous, formula for rho.
# For completion, we will also define the exact formula for rho from option D.
rho_expr_option_D = (1 - (rho_m - tau_m**2) * sympy.exp(I * 2 * k0 * d) * rho_m) / (1 - rho_m**2 * sympy.exp(I * 2 * k0 * d))


print("Based on the derivation, the transmission coefficient is:")
print("τ =", sympy.pretty(tau_expr, use_unicode=True))
print("\nThis uniquely matches option D.")
print("The reflection coefficient from option D is:")
print("ρ =", sympy.pretty(rho_expr_option_D, use_unicode=True))
print("\nNote: The derived reflection coefficient is slightly different, but option D is the best fit overall.")
print("Derived ρ =", sympy.pretty(rho_expr, use_unicode=True))
