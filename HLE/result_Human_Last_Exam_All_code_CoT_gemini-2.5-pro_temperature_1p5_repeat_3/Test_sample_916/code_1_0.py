import sympy
from sympy import Symbol, sin, cos, pi, simplify, Eq, pretty_print

# Define the symbols
f_x = Symbol('f_x(t)')
R = Symbol('R')
N = Symbol('N')
N0 = Symbol('N_0')
I0 = Symbol('I_0')
i0 = Symbol('i_0')
omega = Symbol('omega')
t = Symbol('t')
g = Symbol('g')
mu0 = Symbol('mu_0')
alpha_T = Symbol('alpha_T')
T = Symbol('T')
T0 = Symbol('T_0')
Bs = Symbol('B_s')

# --- Build the expression step-by-step based on the derivation ---

# 1. Time-varying current
i_t = i0 * sin(omega * t)

# 2. Temperature corrected permeability numerator
mu_temp_num = mu0 * (1 - alpha_T * (T - T0))

# 3. Saturation term in the denominator
saturation_term = 1 + (mu0 * N0 * I0) / (g * Bs)

# 4. Numerator of the force expression (excluding the permeability part)
force_numerator_main = -2 * pi * R * N * N0 * I0 * i_t

# 5. Denominator of the force expression
force_denominator = g**2 * saturation_term

# 6. Full force expression by combining numerator and permeability
# f_x = (force_numerator_main * mu_temp_num) / force_denominator
# This is equivalent to Choice B
final_expression = force_numerator_main / g**2 * mu_temp_num / saturation_term

# Create a Sympy Equation object for pretty printing
final_equation = Eq(f_x, final_expression)

# Print the final equation in a readable format
print("The instantaneous force f_x(t) is given by the equation:")
pretty_print(final_equation)

print("\nLet's format the expression to exactly match option B:")
# To match the structure of option B, group terms differently
option_B_expr = -2*pi*R*N * (mu_temp_num * N0 * I0 * i0 * sin(omega*t)) / (g**2 * saturation_term)
option_B_equation = Eq(f_x, option_B_expr)
pretty_print(option_B_equation)