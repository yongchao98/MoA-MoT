import sympy

# Define symbols for the expression
i = sympy.I
pi = sympy.pi
zeta = sympy.functions.special.zeta.zeta

# Define the components of the analytical solution
s = sympy.S(3)/2
q_real = sympy.S(3)/2
q_imag = -sympy.S(2)/3
q = q_real + i * q_imag

# Construct the full expression for the integral I
# I = 2 * i * sqrt(pi) * Re(zeta(s, q))
# where Re(zeta(s, q)) is the real part of the Hurwitz zeta function.
# The following code prints this result.

print("The analytical value of the integral is I, where:")
# The equation for I
equation_I = sympy.Eq(sympy.Symbol('I'), 2 * i * sympy.sqrt(pi) * sympy.Re(zeta(s, q)))

# Print the equation
# We format it to make it more readable, printing each component.
coeff_val = 2
func_name = "Re(zeta)"
param1_val = s
param2_val_real = q_real
param2_val_imag_abs = abs(q_imag)

print(f"I = {coeff_val} * i * sqrt(pi) * {func_name}({param1_val}, {param2_val_real} - {param2_val_imag_abs}*i)")