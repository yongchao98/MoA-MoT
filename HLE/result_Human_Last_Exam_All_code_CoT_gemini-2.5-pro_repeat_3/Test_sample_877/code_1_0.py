import math

# The function h(x) is derived from the equation of the separatrix that
# contains the saddle point (0, 1/2). This separatrix defines the boundary
# for the basin of attraction of the stable equilibrium.
# The equation is a^2 = h(b).
#
# Our derived function is h(x) = 4x^2 + 2x - 2 + 2x*ln(2x).
#
# The problem asks to output each number in the final equation.
# We will print the function's expression with its coefficients.

coeff_x_squared = 4
coeff_x = 2
constant_term = -2
coeff_log_term = 2
coeff_in_log = 2

print("The function h(x) is determined by the equation for the separatrix, which separates bounded and unbounded solutions.")
print("The condition given, -sqrt(h(b(0))) < a(0) < 0, describes the set of initial conditions inside this boundary.")
print("Under the standard interpretation that such a system has implicit dissipation, these initial conditions lead to a(t) -> 0.")
print("\nThe function h(x) is:")
print(f"h(x) = {coeff_x_squared}*x^2 + {coeff_x}*x + ({constant_term}) + {coeff_log_term}*x*ln({coeff_in_log}*x)")
