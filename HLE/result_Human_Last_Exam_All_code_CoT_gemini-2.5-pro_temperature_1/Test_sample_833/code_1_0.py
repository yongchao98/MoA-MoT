import sympy

# Define the symbols
# u0 represents the constant value around which we linearize.
u0 = sympy.Symbol('u0')

# The expression for the lower bound from the heuristic analysis
# is -u0 * (1 - u0) / 2.
expression = -u0 * (1 - u0) / 2

# We need to find the minimum value of this expression for u0 in [0, 1].
# To do this, we can take the derivative with respect to u0 and set it to 0.
derivative = sympy.diff(expression, u0)

# Solve for u0 where the derivative is 0 to find critical points.
critical_points = sympy.solve(derivative, u0)
# The critical point is u0 = 1/2.
u0_min = critical_points[0]

# Substitute this value back into the expression to find the minimum value.
min_value = expression.subs(u0, u0_min)

# The result is the lower bound 'a'.
a = min_value

print(f"The expression for the approximate lower bound is: {expression}")
print(f"The derivative with respect to u0 is: {derivative}")
print(f"The critical point is at u0 = {u0_min}")
print(f"The minimum value (lower bound a) is: {a}")
print("Final Answer Equation:")
# The equation is a <= expression
# We found a = -1/8
print(f"{a} <= (d/dt + (1-2*u)*u_bar*d/dx)*u_bar")
