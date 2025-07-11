import sympy

# Let P be the probability of escaping.
# Based on the analysis, we derived a linear equation for P.
# (5/9) * P = 1/6
# We can represent P as a symbolic variable to show the equation.
P = sympy.Symbol('p')

# The equation is 5/9 * P - 1/6 = 0
lhs = sympy.Rational(5, 9) * P
rhs = sympy.Rational(1, 6)

# We want to present the equation in a simple integer form.
# Multiply by the least common multiple of the denominators (9 and 6), which is 18.
# 18 * (5/9) * P = 18 * (1/6)
# 10 * P = 3
eq_lhs_coeff = 10
eq_rhs_val = 3

# Solve for P
p_val = sympy.solve(eq_lhs_coeff * P - eq_rhs_val, P)[0]

# Print the final equation and the result
print("Let p be the probability that the marble escapes.")
print(f"The derived equation for p is: {eq_lhs_coeff} * p = {eq_rhs_val}")
print(f"The solution for p is {p_val.p}/{p_val.q}")
print(f"Final Answer: {float(p_val)}")