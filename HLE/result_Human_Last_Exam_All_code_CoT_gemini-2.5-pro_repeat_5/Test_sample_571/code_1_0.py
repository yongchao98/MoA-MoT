import sympy

# Define the symbols for the variables
x, y = sympy.symbols('x y')
a, b = -2, -1

# The trace formula we derived
q, z_tr = sympy.symbols('q z_tr')
trace_formula = q**(-1) * z_tr**2 - (1 - q**(-1)) * z_tr

# The substitution rules
q_sub = x**a
z_sub = x**b * y

# Apply the substitution
poly = trace_formula.subs([(q, q_sub), (z_tr, z_sub)])
poly_simplified = sympy.simplify(poly)

# Print the resulting polynomial
# We manually expand to show each number in the final equation
# poly_simplified is y**2 - y*(x**(-1) - x) = 1*y**2 - (1*x**-1 - 1*x)*y
term1_coeff = 1
term2_coeff_x_inv = 1
term2_coeff_x = -1

print(f"The resulting polynomial is P(x,y) = {term1_coeff}*y**2 - ({term2_coeff_x_inv}*x**-1 + ({term2_coeff_x})*x)*y")
print(f"Simplified: P(x,y) = {poly_simplified}")