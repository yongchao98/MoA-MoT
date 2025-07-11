import math

# The function h(x) is derived from the equation for the stable manifold
# of the equilibrium point (a,b) = (0, 1/2).
# Let x = b(0). The derived equation for the manifold is:
# a(0)^2 = 4*b(0)^2 - 6*b(0) + 2 + 2*b(0)*ln(2*b(0))
# Thus, h(x) is the right-hand side of this equation.
# h(x) = 4*x^2 - 6*x + 2 + 2*x*ln(2*x)

# We define the coefficients of the function h(x)
coeff_x_squared = 4
coeff_x = -6
constant_term = 2
coeff_x_log = 2
coeff_in_log = 2

# We format the function as a string to display it clearly,
# including every number as requested.
h_x_expression = (f"h(x) = {coeff_x_squared}*x^2 + ({coeff_x})*x + {constant_term} + "
                  f"{coeff_x_log}*x*ln({coeff_in_log}*x)")

print("The function h(x) that defines the condition for convergence is:")
print(h_x_expression)