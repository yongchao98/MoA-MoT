import sympy

# Define the symbols
n, theta, omega, c, epsilon_0, E_x0_i = sympy.symbols(
    'n, theta, omega, c, epsilon_0, |E_{x0}^i|'
)

# Common denominator term
k_0 = omega / c
alpha_factor = sympy.sqrt(n**2 * sympy.sin(theta)**2 - 1)
common_denominator = 2 * k_0 * (n**2 - 1) * ((n**2 + 1)*sympy.sin(theta)**2 - 1) * alpha_factor

# Numerator for Electric Field Energy
numerator_E = n**2 * (2*n**2 * sympy.sin(theta)**2 - 1)

# Numerator for Magnetic Field Energy
numerator_H = n**2 * (n**2 * sympy.sin(theta)**2 - 1)

# Full expression for Electric Field Energy
Energy_E_expr = (numerator_E / common_denominator) * epsilon_0 * E_x0_i**2

# Full expression for Magnetic Field Energy
Energy_H_expr = (numerator_H / common_denominator) * epsilon_0 * E_x0_i**2

# Print the final equations in a readable format
print("The total time-averaged stored energy per unit area is given by:")
print("-" * 60)
print("Energy in E field = ")
sympy.pprint(Energy_E_expr, use_unicode=True)
print("\n" + "-" * 60)
print("Energy in H field = ")
sympy.pprint(Energy_H_expr, use_unicode=True)
