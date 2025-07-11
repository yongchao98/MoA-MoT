import sympy

# Define variables
x, y, a, b = sympy.symbols('x y a b')
M = sympy.Symbol('M', positive=True) # M is a large positive constant

# Define the two inequalities
# Inequality 1: y >= M*b - M*a - M
# Inequality 2: y >= x - M*b - M*a
ineq1_lhs = y
ineq1_rhs = M*b - M*a - M

ineq2_lhs = y
ineq2_rhs = x - M*b - M*a

# Format the output as requested
print(f"{ineq1_lhs} >= {sympy.expand(ineq1_rhs)}")
print(f"{ineq2_lhs} >= {sympy.expand(ineq2_rhs)}")
