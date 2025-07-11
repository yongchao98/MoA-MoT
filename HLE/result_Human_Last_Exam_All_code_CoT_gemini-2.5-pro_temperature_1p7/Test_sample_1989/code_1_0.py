import sympy as sp

# Define the symbols
r, theta, A, B = sp.symbols('r, theta, A, B', real=True)
C = sp.Symbol('C') # Represents a constant

# The original asymptotic behavior for A=B=0
omega_0_expr = (1/sp.sqrt(r)) * sp.exp(-r*(1-sp.cos(theta)))

# The corrector term due to non-zero A and B
corrector_expr = r**A * sp.exp(B*theta)

# The new asymptotic behavior
omega_expr = corrector_expr * omega_0_expr

# To match the requested format r^(power) * exp(argument)
power_expr = A - sp.Rational(1, 2)
exp_argument_expr = B*theta - r*(1-sp.cos(theta))

final_expr = C * r**power_expr * sp.exp(exp_argument_expr)

# Print the components of the final expression
# The problem asks to output each number in the final equation.
# In the expression C * r**(A - 1/2) * exp(B*theta - r*(1-cos(theta))), the numbers are -1/2, -1, 1.

print("The large-distance behavior of omega(r, theta) is proportional to:")
print(f"r**({power_expr}) * exp({exp_argument_expr})")

# Extract and print the numbers in the final symbolic expression
print("\nLet's break down the expression C * r**(A - 1/2) * exp(B*theta - r*(1-cos(theta))):")

print(f"The power of r is A + ({sp.latex(power_expr.args[1])})")
b_theta_term = exp_argument_expr.args[0]
r_term = exp_argument_expr.args[1]

r_coeff = r_term.args[0]
cos_term = r_term.args[1].args[1]
cos_coeff = cos_term.args[0]

one_minus_cos_term = r_term.args[1]
one_in_term = one_minus_cos_term.args[0]


print(f"The coefficient of theta in the exponent is B.")
print(f"The term proportional to r in the exponent is ({sp.latex(r_coeff)})*r*(...)")
print(f"Inside the parenthesis with r, we have the number {sp.latex(one_in_term)} and the coefficient of cos(theta) which is ({sp.latex(cos_coeff)})")
print(f"So the numbers are: -1/2 (from the power of r), and in the exponent -1 (coefficient of r), 1 (the constant term), and -1 (coefficient of cos(theta)) which is multiplied by -1 again.")
