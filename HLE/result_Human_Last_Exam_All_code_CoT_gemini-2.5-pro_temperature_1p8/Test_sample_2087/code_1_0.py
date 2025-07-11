import sympy as sp

# Define the symbolic variables used in the problem statement
x = sp.Symbol('x')
# Using functions for clarity in the final expression
F_Xi = sp.Function('F_{X_i}')
I_Xi = sp.Function('I_{X_i}')
mu_Xi = sp.Symbol('\\mu_{X_i}')

# Construct the numerator of the expression
numerator = x * F_Xi(x) - I_Xi(x)

# The final expression for the limiting CDF
final_expression = numerator / mu_Xi

# To satisfy the "output each number in the final equation" instruction,
# we will print the final equation in a structured way.
# LHS stands for Left Hand Side, RHS for Right Hand Side.

# Let's represent the left hand side symbolically too.
t = sp.Symbol('t')
X = sp.Function('X')
F_Xt = sp.Function('F_{X(t)}')

lhs = "lim_{t->\infty} " + str(F_Xt(x))
rhs = str(final_expression).replace('F_{X_i}', 'F_{X_i}').replace('I_{X_i}', 'I_{X_i}')

# Manually format for better readability
rhs_formatted = f"({x}*F_{{X_i}}({x}) - I_{{X_i}}({x})) / {mu_Xi}"

# Printing the final equation
print(f"The expression for the limiting CDF is:")
# Python's print will output all the parts of the final constructed string.
print(f"lim_{{t->\infty}} F_{{X(t)}}(x) = {rhs_formatted}")
