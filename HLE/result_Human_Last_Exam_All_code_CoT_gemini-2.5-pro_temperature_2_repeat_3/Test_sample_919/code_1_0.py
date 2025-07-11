import sympy as sp

# Define the symbols
K0, mu0, mu, a, d, y = sp.symbols('K_0 mu_0 mu a d y', real=True, positive=True)
ix = sp.Symbol('i_x')

# The expression for force per unit area from option C
force_per_area_numerator = (mu0 / 2) * K0**2 * sp.sin(a*y)**2
force_per_area_denominator = (sp.cosh(a*d) + (mu0/mu) * sp.sinh(a*d))**2
force_per_area = (force_per_area_numerator / force_per_area_denominator) * ix

# Print the chosen equation in a readable format
print("The force per unit y-z area on the x = d interface is given by:")
# Use sp.pretty_print for a more mathematical layout
sp.pretty_print(sp.Eq(sp.Symbol('f/area'), force_per_area), use_unicode=True)
print("\nThis corresponds to Answer Choice C.")
