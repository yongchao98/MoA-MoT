import sympy

# The problem reduces to finding the minimum value of the expression f(u) = -0.5 * u * (1 - u)
# for u in the interval [0, 1].

# We define the variable u and the function f(u)
u = sympy.symbols('u')
f_u = -0.5 * u * (1 - u)

# Find the derivative of f(u) with respect to u
f_prime_u = sympy.diff(f_u, u)

# Solve for the critical points by setting the derivative to zero
critical_points = sympy.solve(f_prime_u, u)

# The critical point is the value of u that minimizes the expression
u_min = critical_points[0]

# Calculate the minimum value of the function, which is the lower bound 'a'
a = f_u.subs(u, u_min)

# Print the final equation with the values substituted
c1 = 0.5
c2 = 1.0
u_val = float(u_min)
print(f"The lower bound 'a' is calculated from the minimum of the function f(u) = -c1 * u * (c2 - u).")
print(f"The minimum occurs at u = {u_val}.")
print(f"Final Equation: a = -{c1} * {u_val} * ({c2} - {u_val})")
print(f"a = {float(a)}")