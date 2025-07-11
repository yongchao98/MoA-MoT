import sympy

# The condition for the integral to diverge is that the exponent of the integrand
# in the radial direction must be greater than or equal to -1.
# The exponent is k = power_of_r_from_measure - p * decay_exponent_of_I.
#
# On the (a_1, a_3) subspace, the decay exponent of I is 1/2 (with respect to a_3).
# The measure contributes a power of 1.
# So the exponent k is 1 - p * (1/2).
k, p = sympy.symbols('k p')

# Set up the equation for the critical exponent
power_of_r = 1
decay_exponent = 1/2
k_expr = power_of_r - p * decay_exponent

# The divergence occurs when k >= -1. The critical point is k = -1.
critical_exponent_val = -1
equation = sympy.Eq(k_expr, critical_exponent_val)

# Solve for p
solution = sympy.solve(equation, p)
p_critical = solution[0]

# Print out the reasoning and the final equation
print("The condition for the integral to diverge is when the exponent k is greater than or equal to -1.")
print(f"The exponent k is given by the power from the measure ({power_of_r}) minus p times the decay rate ({decay_exponent}).")
print(f"So, k = {power_of_r} - p * {decay_exponent}")
print(f"The critical case for divergence is k = {critical_exponent_val}.")
print("\nThe final equation to solve for p is:")
print(f"{equation.lhs.args[0]} + {equation.lhs.args[1].args[0]}*p*{equation.lhs.args[1].args[1]} = {equation.rhs}")

print(f"\nSolving for p gives the largest value for which the function is not in L^p.")
print(f"The largest value of p is: {p_critical}")
