import sympy

# Since no numerical values are provided, we will use symbolic variables to represent the formula.
# We will define all the variables from the problem description.
R, N, N0, I0, i0, omega, t, g, mu0, alpha_T, T, T0, Bs = sympy.symbols(
    'R N N_0 I_0 i_0 omega t g mu_0 alpha_T T T_0 B_s'
)

# The instantaneous force f_x(t) from choice B is constructed as follows:
# Numerator part
numerator = -2 * sympy.pi * R * N * mu0 * (1 - alpha_T * (T - T0)) * N0 * I0 * i0 * sympy.sin(omega * t)

# Denominator part
denominator = g**2 * (1 + (mu0 * N0 * I0) / (g * Bs))

# The complete expression for the force
force_expression = numerator / denominator

# To present the final formula clearly, we print it using sympy's pretty print.
# This fulfills the requirement to "output each number in the final equation" by showing the symbolic form.
print("f_x(t) = ", end="")
sympy.pprint(force_expression, use_unicode=True)
